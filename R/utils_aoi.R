#' @title internal functions to work with tree list data within a specified area of interest (AOI) or "domain" in the terminology of really smart fire modellers.
#'
#' @description
#' internal functions to work with tree list data within an AOI.
#' For example, to use the tree list (e.g. as exported by [raster2trees()]) within the QUIC-Fire modelling tool some additional columns are required
#' and one needs to define a study area and align the DTM and tree list with this domain.
#'
#' @param tree_list data.frame. A data frame with the columns `treeID`, `tree_x`, and `tree_y`.
#' If an `sf` class object with POINT geometry (see [sf::st_geometry_type()]), the program will use the data "as-is" and only require the `treeID` column.
#' Other required columns include:
#' * `crown_area_m2`, `tree_height_m` (e.g. as exported by [raster2trees()])
#' * `tree_cbh_m` (e.g. as exported by [trees_cbh()])
#' * and one of `dbh_cm`, `dbh_m`, or  `basal_area_m2` (e.g. as exported by [trees_dbh()])
#' @param crs string. A crs string as returned from [sf::st_crs()] or the EPSG code of the x,y coordinates.
#' Defaults to the crs of the `tree_list` data if of class "sf".
#' @param study_boundary sf. The boundary of the study area to define the area of interest which may extend beyond the space with trees.
#' This must be an sf class object with a single record. If you need to get the trees within multiple different AOI's, then `purrr::map()` this function.
#' @param bbox_of_aoi logical. Should the study_boundary be transformed to a bounding box instead of it's original shape for determining the trees within the boundary?
#' If set to true, the bounding box is created prior to applying the buffer.
#' @param buffer numeric. Buffer to be applied to the study area prior to determining trees within the boundary. 
#' Units are determined by the horizontal CRS settings of the tree_list data or the CRS of the reproject_epsg.
#' @param reproject_epsg numeric. The EPSG code to reproject the data in prior to buffering and clipping.
#' Will determine the projection of the output data.
#'
#' @keywords internal
#'
clip_tree_list_aoi <- function(
  tree_list
  , crs
  , study_boundary
  , bbox_of_aoi = F
  , buffer = 0
  , reproject_epsg = NULL
) {
  # function to clip a tree list to an study_boundary

  ###############################
  # checks
  ###############################
  if(!inherits(tree_list,"sf")){stop("tree_list must be sf class object")}
  if(is.na(sf::st_crs(tree_list))){stop("tree_list does not have a CRS")}
  
  # set reproject_epsg if null
  if(is.null(reproject_epsg)){
    reproject_epsg <- sf::st_crs(tree_list)
  }
  # custom boundary
  study_boundary <- get_custom_aoi(
      study_boundary = study_boundary
      , bbox_of_aoi = dplyr::coalesce(bbox_of_aoi,FALSE)
      , buffer = buffer
      , reproject_epsg = reproject_epsg
    )

  ###############################
  # data transformations
  ###############################
  # convert tree list to spatial points data
  tree_tops <- check_spatial_points(tree_list, crs)
  
  # reproj...if is.null(reproject_epsg) initially, this line has no impact
  tree_tops <- tree_tops %>% sf::st_transform(crs = sf::st_crs(study_boundary))
  
  ###############################
  # intersect
  ###############################
  # intersect based on points but filter original tree list
  tree_list <- tree_list %>% 
    dplyr::slice(
      sf::st_intersects(tree_tops, study_boundary, sparse = F) %>% which()
    )
  
  # return will include the tree list and the spatial aoi used to filter the list
  ret <- list(
      tree_list = tree_list
      , aoi = study_boundary # will be the bbox as an sf if bbox_of_aoi
    )

  if(nrow(tree_list)==0){
    warning("no tree_list trees found within study_boundary bounds")
    return(ret)
  }else{
    return(ret)
  }

}
#######################################################
# intermediate function 2
# get an aoi with/without a buffer and/or as/not as a bbox
#######################################################
get_custom_aoi <- function(
  study_boundary
  , bbox_of_aoi = F
  , buffer = 0
  , reproject_epsg = NULL
) {
  # bounds check
  if(!inherits(study_boundary,"sf")){stop("study_boundary must be sf class object")}
  if(is.na(sf::st_crs(study_boundary))){stop("study_boundary does not have a CRS")}
  if(nrow(study_boundary)!=1){
    stop("study_boundary must only have a single record geometry")
  }
  if(
    !all( sf::st_is(study_boundary, c("POLYGON","MULTIPOLYGON")) )
  ){
    stop("study_boundary must contain POLYGON type geometry only")
  }
  
  # check epsg
  if(inherits(reproject_epsg,"crs")){
    # transform
    study_boundary <- study_boundary %>% sf::st_transform(reproject_epsg)
  }else if(
    !is.null(reproject_epsg)
  ){
    # if character
    if(is.character(reproject_epsg)){
      reproject_epsg <- readr::parse_number(reproject_epsg)
    } # if numeric, just use it

    # set the crs to use 
    this_crs <- sf::st_crs(reproject_epsg)
    if(is.na(this_crs)){
      stop("could not find CRS information for given reproject_epsg. make sure this is just the numeric EPSG code.")
    }

    # transform
    study_boundary <- study_boundary %>% sf::st_transform(reproject_epsg)
  }

  # bbox
  if(bbox_of_aoi){
    study_boundary <- sf::st_bbox(study_boundary) %>% 
      sf::st_as_sfc()
  }
  
  # buff
  # check buff
  if(is.character(buffer)){
    buffer <- readr::parse_number(buffer) 
  }
  if(dplyr::coalesce(buffer,0)>0){
    study_boundary <- sf::st_buffer(study_boundary, buffer)
  }

  return(study_boundary)
}

#######################################################
# intermediate function 3
# function to make the simulation bounds/domain for QUIC-Fire
#######################################################
quicfire_define_domain <- function(
  sf_data
  , horizontal_resolution = 2
  , outdir = tempdir()
){
  # checks
  if(!inherits(sf_data,"sf")){stop("sf_data must be sf class object")}
  if(is.na(sf::st_crs(sf_data))){stop("sf_data does not have a CRS")}

  horizontal_resolution <- as.numeric(horizontal_resolution)
  if(
    is.null(horizontal_resolution) 
    || is.na(horizontal_resolution)
    || !inherits(horizontal_resolution,"numeric")
  ){
    stop("horizontal_resolution must be numeric")
  }

  # outdir
  outdir <- file.path(outdir)
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }

  # Calculate initial bounding box
  bbox <- sf::st_bbox(sf_data)

  # Update bounding box to be rounded to be even (so that it grids nicely for modeling later and so we can get nx and ny)
  # get dimensions
  xmax <- as.numeric(bbox$xmax)
  xmin <- as.numeric(bbox$xmin)
  ymin <- as.numeric(bbox$ymin)
  ymax <- as.numeric(bbox$ymax)
  # get length and width
  width <- xmax - xmin
  length <- ymax - ymin
  # Ensure width is divisible by horizontal resolution
  if (width %% horizontal_resolution != 0) {
    xmax <- horizontal_resolution * round(xmax/horizontal_resolution)
    xmin <- horizontal_resolution * round(xmin/horizontal_resolution) - horizontal_resolution
    }
  # Ensure length is divisible by horizontal resolution
  if (length %% horizontal_resolution != 0) {
    ymax <- horizontal_resolution * round(ymax/horizontal_resolution)
    ymin <- horizontal_resolution * round(ymin/horizontal_resolution) - horizontal_resolution
  }
  # recalculate width and length (m)
  width <- xmax - xmin
  length <- ymax - ymin

  #calc domain width and length
  nx <- width/horizontal_resolution #horizontal
  ny <- length/horizontal_resolution #horizontal

  #new bbox rounded up and down and write ff_domain.geojson in WGS84 projection for GUI
  bbox <- sf::st_bbox(
      c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)
      , crs = sf::st_crs(sf_data)
    ) %>% 
    sf::st_as_sfc() %>% 
    sf::st_transform(crs=4326)

  # write it
  # writes to geojson in WGS84
  sf::st_write(
    bbox
    , file.path(outdir,"Lidar_Bounds.geojson")
    , driver = "GeoJSON"
    , delete_dsn = TRUE
  )

  message(paste0(
    "exported QUIC-Fire domain to:"
    , "\n ........ "
    , file.path(outdir,"Lidar_Bounds.geojson")
  ))
  
  return(list("xmin" = xmin, "ymin" = ymin,"nx" = nx, "ny"= ny, "width" = width, "length" = length, "bbox"=bbox))

}

#######################################################
# intermediate function 6
# function to make a Fortran topography file for QUIC-Fire O_o
#######################################################
quicfire_dtm_topofile <- function(
  dtm_rast
  , horizontal_resolution = 2
  , study_boundary = NULL
){
  # open dtm
  if(inherits(dtm_rast, "SpatRaster")) {
    dtm_rast <- dtm_rast
  }else if(inherits(dtm_rast,"raster")){
    dtm_rast <- terra::rast(dtm_rast)
  }else if(
    inherits(dtm_rast, "character")
  ){
    if(
      !any(stringr::str_ends(dtm_rast, ".*\\.(tif|tiff)$")) 
      && dir.exists(file.path(dtm_rast))
    ){
      # try to read directory for dtm_rast files
      fls <- list.files(
        file.path(dtm_rast)
        , pattern = "dtm_.*\\.tif$"
        , full.names = TRUE
      )
      # stop it if no files
      if(length(fls)<1){
        stop(paste0(
          "no DTM raster file found at: "
          , "\n .... "
          , file.path(dtm_rast)
        ))
      }
      # read it
      dtm_rast <- terra::rast(fls[1]) ## only reads one file...but what if multiple?

    }else if(any(stringr::str_ends(dtm_rast, ".*\\.(tif|tiff)$"))){
      # read it
      dtm_rast <- stringr::str_subset(dtm_rast, pattern = ".*\\.(tif|tiff)$")[1] %>%
        ## only reads one file...but what if multiple?
        terra::rast()
    }else{
      stop("this is not a readabile dtm_rast")
    }
  }else{
    stop("dtm_rast is not of class terra")
  }

  #check the crs
  if(is.na(terra::crs(dtm_rast))){
    stop("the crs for dtm_rast is NA")
  }
  
  # we need to make the raster in the desired resolution
  horizontal_resolution <- as.numeric(horizontal_resolution)
  if(
    is.null(horizontal_resolution) 
    || is.na(horizontal_resolution)
    || !inherits(horizontal_resolution,"numeric")
  ){
    stop("horizontal_resolution must be numeric")
  }
  dtm_rast <- adjust_raster_resolution(dtm_rast, target_resolution = horizontal_resolution)

  #clip dtm to domain bounds
  Lidar_Bounds <- st_transform(st_read(paste0(project_path,"Lidar_Bounds.geojson")), crs = epsg_data) #read shapefile we made earlier
  clipped_dtm <- crop(dtm, Lidar_Bounds)
  #sum(is.na(values(clipped_dtm)))#count NANS
  
  #write the clipped dtm to a tif file just to be nice
  writeRaster(clipped_dtm, paste0(project_path,"dtm_Clipped"), format = "GTiff",overwrite=TRUE) #write clipped dtm raster to tif to check


  ### Write to FORTRAN File in format needed for QUIC-Fire (and FIRETEC)

  # Flip the raster over the y-axis
  raster_data <- flip(clipped_dtm, direction = 'y')
  
  # Extract the raster values
  values <- getValues(raster_data)
  
  # Open a connection to the unformatted Fortran file
  fortran_file <- file(paste0(project_path,"topo.dat"), "wb")
  
  # Write the data to the file
  writeBin(charToRaw("BRUH"), fortran_file) #Topo requires a header of 4 bits... this works and brings me joy
  writeBin(values, fortran_file, size = 4)  # Assuming 4-byte (32-bit) floats
  
  # Close the file connection
  close(fortran_file)
  print('Written topo.dat to fortran file')
}
