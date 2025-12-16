#' @title tools for working with projections, crs, epsg, and whatnot
#'
#' @param x sf. or any object of a class that can be read by sf::st_crs such as "LAS" class
#'
#' @return A data.frame with the units and scale factor for the vertical projection
#'
#' @keywords internal
#'
get_vertical_crs <- function(x) {
  xcrs <- sf::st_crs(x)
  if (is.na(xcrs)){
    # stop("No CRS defined...try setting the parameter `new_crs` if known")
    return(NA)
  }

  wkt <- sf::st_as_text(xcrs)

  if (!grepl("COMPD_CS", wkt)) {
    # Should just be a horizontal CRS - nothing
    message("no vertical crs found")
    return(NA)
  }

  # extract the VERT_CS component
  if(!stringr::str_detect(wkt,"VERT_CS")){
    message("no vertical crs found")
    return(NA)
  }

  # get the vertical component
  # might also be VERT_DATUM ???
  i <- regexpr("VERT_CS\\[", wkt)
  wkt <- base::substring(wkt, i)

  # find the units
  if(!stringr::str_detect(wkt,"UNIT")){
    message("no vertical crs found")
    return(NA)
  }

  i <- regexpr("UNIT\\[", wkt)
  wkt <- base::substring(wkt, i) %>% stringr::str_remove("UNIT\\[")

  # work with the unit part
  unit_part <- stringr::str_extract(wkt, "[^\\]]+")
  unit_df <-
    read.csv(text = unit_part, header = F, colClasses = "character") %>%
    dplyr::mutate(
      units = stringr::str_extract(
        tolower(V1)
        , paste(c("foot","feet","meter","metere","metre"),collapse = "|")
      )
      , scale_factor = readr::parse_number(V2)
    ) %>%
    dplyr::select(units,scale_factor) %>%
    dplyr::slice(1)

  if(nrow(unit_df)==1){
    return(unit_df)
  }else{
    return(NA)
  }

}

###___________________________________________###
# Retrieve the horizontal component of a compound CRS.
# The object x can be an 'sf' package 'crs' object or any
# spatial object from which a CRS can be queried using the
# sf::st_crs function.
###___________________________________________###
get_horizontal_crs <- function(x) {
  xcrs <- sf::st_crs(x)
  if (is.na(xcrs)){
    # stop("No CRS defined...try setting the parameter `new_crs` if known")
    return(NA)
  }

  wkt <- sf::st_as_text(xcrs)

  if (!grepl("COMPD_CS", wkt)) {
    # Should just be a horizontal CRS - simply return it
    xcrs
  } else {
    # Extract the horizontal component
    i <- regexpr("PROJCS\\[", wkt)
    wkt <- base::substring(wkt, i)

    # Match square brackets to discard any trailing
    # component (e.g. the vertical CRS)
    wkt_chars <- base::strsplit(wkt, "")[[1]]
    level <- 1
    k <- base::match("[", wkt_chars)
    while (level > 0) {
      k <- k + 1
      if (wkt_chars[k] == '[') {
        level <- level + 1
      } else if (wkt_chars[k] == ']') {
        level <- level - 1
      }
    }

    wkt <- base::substring(wkt, 1, k)
    # return
    return(sf::st_crs(wkt))
  }
}

###___________________________________________###
# transform raw las x,y points and return a data.frame
# of x,y in new projection which can be combined with
# a Z column for full definition of LAS
# https://gis.stackexchange.com/questions/375735/clean-way-to-convert-from-us-feet-to-meters-with-lidr-lidar-point-cloud
# https://github.com/r-lidar/lidR/issues/372
###___________________________________________###
st_transform_las_xy <- function(las, new_epsg_code) {
  stopifnot(inherits(las, "LAS"))
  # get epsg number
  e <- check_epsg_code(new_epsg_code = new_epsg_code)
  # do the work and return data.frame with X,Y that can be converted to las after adding Z
  if(lidR::is.empty(las)){return(NULL)}
  return_df <- las@data %>%
    dplyr::select(X,Y) %>%
    sf::st_as_sf(
      coords = c("X", "Y")
      , crs = lidR::st_crs(las)
      , remove = F
    ) %>%
    dplyr::rename(orig_x = X, orig_y = Y) %>%
    sf::st_transform(crs = e) %>%
    dplyr::mutate(
      X = sf::st_coordinates(.)[,1]
      , Y = sf::st_coordinates(.)[,2]
    ) %>%
    dplyr::select(X,Y) %>%
    sf::st_drop_geometry()
  # return
  return(return_df)
}

###___________________________________________###
# transform raw las z points and return a data.frame
# of z scaled which can be combined with
# a x,y column for full definition of LAS
# https://gis.stackexchange.com/questions/375735/clean-way-to-convert-from-us-feet-to-meters-with-lidr-lidar-point-cloud
# https://github.com/r-lidar/lidR/issues/372
###___________________________________________###
st_transform_las_z <- function(
  las
  , scale_factor # this gets multiplied by current z
){
  stopifnot(inherits(las, "LAS"))
  # get epsg number
  if(inherits(scale_factor,"character")){
    f <- readr::parse_number(scale_factor)
  }else{
    f <- scale_factor
  }
  if(
    is.na(f) || is.null(f) ||
    identical(f,character(0)) ||
    !inherits(f,"numeric")
  ){
    stop("could not detect numeric scale factor")
  }

  # do the work and return data.frame with z that can be converted to las after adding x,y
  if(lidR::is.empty(las)){return(NULL)}
  return_df <- las@data %>%
    dplyr::select(Z) %>%
    dplyr::mutate(
      Z = Z*f
    )
  # return
  return(return_df)
}

###___________________________________________###
# combine xy and z data frames to make a new las
# uses output from from st_transform_las_xy()
#  and from st_transform_las_z()
# make sure to set new_epsg_code same as used in st_transform_las_xy()
###___________________________________________###
combine_xy_z_make_las <- function(
  xy_df # direct from st_transform_las_xy() without any transformation
  , z_df # direct from st_transform_las_z() without any transformation
  , new_epsg_code # same as used in st_transform_las_xy()
  , file # A character string naming an output file lidR::writeLAS()
  , index = FALSE # boolean. Also write a lax file to index the points in the files lidR::writeLAS()
) {
  stopifnot(inherits(xy_df, "data.frame"))
  stopifnot(inherits(z_df, "data.frame"))
  ## check for cols
  check_df_cols_all_missing(
    xy_df
    , col_names = c("X","Y")
    , all_numeric = T
  )
  check_df_cols_all_missing(
    z_df
    , col_names = c("Z")
    , all_numeric = T
  )
  # get epsg number
  e <- check_epsg_code(new_epsg_code = new_epsg_code)

  # combine xyz for small sample just to set up header
  new_las <- dplyr::bind_cols(
      xy_df %>% dplyr::select(X,Y) %>% dplyr::slice_head(n=11)
      , z_df %>% dplyr::select(Z) %>% dplyr::slice_head(n=11)
    ) %>%
    lidR::LAS(crs = sf::st_crs(e))

  #### need a header
  # lidR::header(new_las)
  # rlas::header_create()
  # lidR::LASheader()
  las_header <- lidR::header(new_las) # lidR::LASheader()
  las_header@VLR <- list() # Erase VLR previously written
  las_header@PHB[["Global Encoding"]][["WKT"]] <- TRUE
  las_header@PHB[["Version Minor"]] <- 4L
  las_header@PHB[["Header Size"]] <- 375L
  las_header@PHB[["Offset to point data"]] <- 375L
  lidR::crs(las_header) <- paste0("EPSG:", e) # not sure if this is needed since we set wkt below
  lidR::wkt(las_header) <- sf::st_crs(e)$wkt
  # las_header

  # combine xyz for full set with header
  new_las <- dplyr::bind_cols(
      xy_df %>% dplyr::select(X,Y)
      , z_df %>% dplyr::select(Z)
    ) %>%
    lidR::LAS(crs = sf::st_crs(e), header = las_header)

  ### warning: Detection of quantization errors for X
  ### see: las_tools.R: https://github.com/r-lidar/lidR/blob/0a09bcb898ce58634e93b200ea0c6db5345b8ae8/R/las_tools.R
  ### Rescale and reoffset recompute the coordinates with
  ### new scales and offsets according to LAS specification
  # new_las <- lidR::las_rescale(new_las)
  # new_las <- lidR::las_reoffset(new_las)
  ### !!!! including lidR::las_rescale & lidR::las_reoffset did not remove the warning...
  ### !!!! ...mapping of las looks good so let's roll with it even if warnings
  ### !!!! based on: https://github.com/r-lidar/lidR/blob/0a09bcb898ce58634e93b200ea0c6db5345b8ae8/R/methods-LAS.R
  ### !!!! lidR::LAS() defaults the x,y offset based on the minimum in the data set and the z offset to "0"
  ### !!!! lidR::LAS() defaults the scale to 0.001 (1 mm if projection is m) for x,y,z
  ### !!!! based on these defaults and an x value of 900567.58294561 (m) which is also the minimum value, x will be stored as:
  ### !!!! (x - offset)/scale = floor( (900567.58294561 - floor(900567.58294561))/0.001 ) = 582 (the decimal part to the mm)
  ### !!!! with an offset set to floor(900567.58294561) = 900567
  ### !!!! ...this all just helps to save on memory as all values are stored as integers and is good enough for our purposes

  # write
  if(!missing(file) && inherits(file,"character")){
    lidR::writeLAS(new_las, file = file, index = index)
  }

  # return las
  return(new_las)
}

###___________________________________________###
# combine all steps above to make a new transformed las
###___________________________________________###
st_transform_las <- function(
    las
    , new_epsg_code
    , file # A character string naming an output file lidR::writeLAS()
    , index = FALSE # boolean. Also write a lax file to index the points in the files lidR::writeLAS()
) {
  # could move to parameters
  transform_only_z_feet <- T

  # checks
  stopifnot(inherits(las, "LAS"))
  if(lidR::is.empty(las)){return(NULL)}
  # get epsg number
  e <- check_epsg_code(new_epsg_code = new_epsg_code)

  # transform xy
  xy_df <- st_transform_las_xy(las, new_epsg_code = e)
  if(dplyr::coalesce(nrow(xy_df),0)==0){return(NULL)}

  # get vertical crs
  vert_crs_df <- get_vertical_crs(las)

  # transform z (or don't)
  if(
    inherits(vert_crs_df, "data.frame") && # df and not NA
    dplyr::coalesce(nrow(vert_crs_df),0)==1 && # has records
    !is.na(vert_crs_df$scale_factor[1]) && # has scale_factor
    !identical(vert_crs_df$scale_factor[1],1) && # scale_factor!=1
    ( # is feet
      transform_only_z_feet &
      !is.na(vert_crs_df$units[1]) &
      stringr::str_detect(
        vert_crs_df$units[1]
        ,paste(c("foot","feet"),collapse = "|")
      )
    )
  ){
    z_df <- st_transform_las_z(las, scale_factor = vert_crs_df$scale_factor[1])
  }else{ # no transform
    z_df <- las@data %>% dplyr::select(Z)
  }

  # combine
  new_las <- combine_xy_z_make_las(
    xy_df = xy_df # direct from st_transform_las_xy() without any transformation
    , z_df = z_df # direct from st_transform_las_z() without any transformation
    , new_epsg_code = e # same as used in st_transform_las_xy()
    , file = file # A character string naming an output file lidR::writeLAS()
    , index = index # boolean. Also write a lax file to index the points in the files lidR::writeLAS()
  )
  if(lidR::is.empty(new_las)){return(NULL)}
  # return
  return(new_las)
}

###___________________________________________###
# st_transform_las on a LASCatalog
# via lidR::catalog_apply()
###___________________________________________###
ctg_st_transform_las <- function(
  chunk
  , new_epsg_code
){
  las <- lidR::readLAS(chunk)
  if(lidR::is.empty(las)){return(NULL)}
  # get epsg number
  e <- check_epsg_code(new_epsg_code = new_epsg_code)
  # transform chunk
  # safe_st_transform_las <- purrr::safely(st_transform_las)
  new_las <- st_transform_las(las = las, new_epsg_code = e)
  # just get the result
  # new_las <- new_las$result
  # return
  return(new_las)
}

###___________________________________________###
# reproject ctg or las with
# st_transform_las()
# always returns a LAScatalog or fails
###___________________________________________###
apply_st_transform_las <- function(
  las # LAScatalog, LAS, dir with las, or las fname
  , outfolder
  , new_epsg_code
) {
  # dir
  if(!dir.exists(outfolder)){
    dir.create(file.path(outfolder), showWarnings = F)
  }

  # check
  new_ctg <- check_las_data(las)
  # set the lascatalog options
  if(inherits(new_ctg, "LAScatalog")){
    # options
    lidR::opt_progress(new_ctg) <- F
    lidR::opt_select(new_ctg) <- "xyz" # 0 enables all extra bytes to be loaded...possibly treeID
    lidR::opt_output_files(new_ctg) <- paste0(file.path(outfolder), "/{*}_reproj")
    # run it
    output_temp <- lidR::catalog_apply(
      ctg = new_ctg
      , FUN = ctg_st_transform_las
      , .options = list(automerge = TRUE)
      # st_transform_las options
      , new_epsg_code = new_epsg_code
    )
    # return ctg
    if(inherits(output_temp, "LAScatalog")){
      output_temp <- lidR::readLAScatalog(output_temp$filename)
    }
  }else if(inherits(new_ctg, "LAS")){
    ofile <- file.path(outfolder, "las_reproj.las")
    # run it
    output_temp <- st_transform_las(las = new_ctg, new_epsg_code = new_epsg_code, file = ofile, index = T)
    # return ctg
    if(
      inherits(output_temp, "LAS") &&
      file.exists(ofile)
    ){
      output_temp <- lidR::readLAScatalog(ofile)
    }
  }else{
    stop("could not reproject")
  }

  if(!inherits(output_temp, "LAScatalog")){
    stop("could not reproject")
  }

  return(output_temp)

}


###___________________________________________###
# adjust the resolution of a raster to be in exactly the target resolution
###___________________________________________###
adjust_raster_resolution <- function(
  raster_object
  , target_resolution
  , fun = mean
  , resample_method = "bilinear"
  , ofile = NULL
) {
  # check if the input is a spatraster object
  if (!inherits(raster_object, "SpatRaster")) {
    stop("Input must be a SpatRaster object.")
  }

  current_resolution <- terra::res(raster_object)[1] # get current resolution (assuming square pixels)
  result_raster <- NULL

  # aggregating (decreasing resolution)
  if (target_resolution > current_resolution) {
    # calculate the aggregation factor
    # we aim for an integer factor for aggregate, but then refine with resample
    fact <- max(1, floor(target_resolution / current_resolution))

    # aggregate the raster
    if(inherits(ofile,"character")){
      aggregated_raster <- terra::aggregate(raster_object, fact = fact, fun = fun, na.rm = TRUE, filename = ofile, overwrite = TRUE)
    }else{
      aggregated_raster <- terra::aggregate(raster_object, fact = fact, fun = fun, na.rm = TRUE)
    }
    result_raster <- aggregated_raster
  }else if(target_resolution < current_resolution) {
    # disaggregating (increasing resolution)

    # calculate the disaggregation factor
    # we aim for an integer factor for disaggregate, but then refine with resample
    fact <- max(1, floor(current_resolution / target_resolution)) # round down to ensure disagg factor is not too large

    # disaggregate the raster
    if(inherits(ofile,"character")){
      disaggregated_raster <- terra::disagg(raster_object, fact = fact, filename = ofile, overwrite = TRUE)
    }else{
      disaggregated_raster <- terra::disagg(raster_object, fact = fact)
    }

    result_raster <- disaggregated_raster

  }else if(target_resolution == current_resolution){
    return(raster_object)
  } else {
    stop("this resolution is unresovable D: ")
  }

  # check if the resulting resolution is exactly the target resolution
  if (abs(terra::res(result_raster)[1] - target_resolution) > 0.0001) { # Using a small tolerance for comparison
    message("the initial aggregation/disaggregation did not result in the exact target resolution. resampling to achieve the precise target resolution.")

    # create a dummy raster with the desired resolution and extent for resampling
    template_raster <- terra::rast(result_raster)
    terra::res(template_raster) <- target_resolution

    # resample the result_raster to the exact target resolution
    if(inherits(ofile,"character")){
      result_raster <- terra::resample(result_raster, template_raster, method = resample_method, filename = tempfile(fileext = ".tif"), overwrite = TRUE)
      terra::writeRaster(result_raster, filename = ofile, overwrite=T)
    }else{
      result_raster <- terra::resample(result_raster, template_raster, method = resample_method)
    }
  }

  return(result_raster)
}

###___________________________________________###
# check epsg code is numeric
###___________________________________________###
check_epsg_code <- function(new_epsg_code) {
  if(inherits(new_epsg_code,"character")){
    e <- readr::parse_number(new_epsg_code)
  }else{
    e <- new_epsg_code
  }
  if(
    is.na(e) || is.null(e) ||
    identical(e,character(0)) ||
    !inherits(e,"numeric")
  ){
    stop("could not detect numeric epsg code")
  }
  return(e)
}

###___________________________________________###
# check_horizontal_crs_is_feet
###___________________________________________###
check_horizontal_crs_is_feet <- function(las) {
  data_crs_params <- las@data %>% sf::st_crs(parameters=T) # wow whaaaaa??? the goods?
  las_crs_params <- get_horizontal_crs(las) %>% sf::st_crs(parameters=T) # wow whaaaaa??? the goods?
  # any matches
  m <- any(
    data_crs_params$ud_unit %>% deparse() %>% stringr::str_detect(paste(c("foot","feet"),collapse = "|")) %>% any()
    , data_crs_params$units_gdal %>% stringr::str_detect(paste(c("foot","feet"),collapse = "|")) %>% any()
    , data_crs_params$proj4string %>% stringr::str_detect("\\+units=us-ft") %>% any()
    # , las_crs_params$wkt %>% stringr::str_detect(paste(c("foot","feet"),collapse = "|"))
    # , las_crs_params$input %>% stringr::str_detect(paste(c("foot","feet"),collapse = "|"))
    , las_crs_params$ud_unit %>% deparse() %>% stringr::str_detect(paste(c("foot","feet"),collapse = "|")) %>% any()
    , las_crs_params$units_gdal %>% stringr::str_detect(paste(c("foot","feet"),collapse = "|")) %>% any()
    , las_crs_params$proj4string %>% stringr::str_detect("\\+units=us-ft") %>% any()
    , na.rm = T
  )
  if(is.na(m) || is.null(m) || identical(m,character(0))){m <- F}
  return(m)
}

###___________________________________________###
# calculate diameter of single polygon
###___________________________________________###
# function to calculate the diamater of an sf polygon that is potentially irregularly shaped
# using the distance between the farthest points
st_calculate_diameter_polygon <- function(polygon) {
  # get the convex hull
  ch <- sf::st_convex_hull(polygon)

  # cast to multipoint then point to get individual vertices
  ch_points <- sf::st_cast(ch, 'MULTIPOINT') %>% sf::st_cast('POINT')

  # calculate the distances between all pairs of points
  distances <- sf::st_distance(ch_points)

  # find the maximum distance, which is the diameter
  diameter <- as.numeric(max(distances,na.rm=T))
  return(diameter)
}
# apply st to sf data
st_calculate_diameter <- function(sf_data) {
  if(!inherits(sf_data,"sf")){stop("st_calculate_diameter() requires polygon sf data")}
  if(
    !all( sf::st_is(sf_data, c("POLYGON","MULTIPOLYGON")) )
  ){
    stop("st_calculate_diameter() requires polygon sf data")
  }

  # get the geometry column name
  geom_col_name <- attr(sf_data, "sf_column")

  # calculate diameter
  # !!rlang::sym() unquotes the geometry column
  return_dta <- sf_data %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(diameter_m = st_calculate_diameter_polygon( !!rlang::sym(geom_col_name) )) %>%
    dplyr::ungroup()
  return(return_dta)
}

###___________________________________________###
# THIS IS OLD AND MAY NOT WORK
# Function to reproject las data using `lidR` but be careful!
# This is inefficient and potentially causes inaccuracies due to transformations (see reference).
# https://gis.stackexchange.com/questions/371566/can-i-re-project-an-las-file-in-lidr
###___________________________________________###
reproject_las <- function(filepath, new_crs = NA, old_crs = NA, outdir = getwd()) {
    if(is.null(new_crs) | is.na(new_crs)){stop("the new_crs must be provided")}
    # read individual file
    las <- lidR::readLAS(filepath)
    if(
      (is.null(sf::st_crs(las)$epsg) | is.na(sf::st_crs(las)$epsg))
      & ( is.null(old_crs) | is.na(old_crs) | old_crs=="" )
    ){
      stop("the raw las file has missing CRS and cannot be transformed. try setting old_crs if known")
    }else{
      # transform if know old crs
      if(
        (is.null(sf::st_crs(las)$epsg) | is.na(sf::st_crs(las)$epsg))
        & !is.null(old_crs) & !is.na(old_crs) & length(as.character(old_crs))>0
      ){
        sf::st_crs(las) <- paste0("EPSG:", old_crs)
      }
      # get filename
      fnm <- filepath %>% basename() %>% stringr::str_remove_all("\\.(laz|las)$")
      new_fnm <- paste0(normalizePath(outdir),"/",fnm,"_epsg",new_crs,".las")
      # reproject
      las <- sf::st_transform(las, paste0("EPSG:", new_crs))
      # write
      lidR::writeLAS(las, file = new_fnm)
      return(new_fnm)
    }
  }
