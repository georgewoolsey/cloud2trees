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
