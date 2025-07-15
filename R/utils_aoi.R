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
  if(
    !inherits(study_boundary,"sf")
    && !inherits(study_boundary,"sfc")
  ){stop("study_boundary must be sf class object")}
  if(is.na(sf::st_crs(study_boundary))){stop("study_boundary does not have a CRS")}
  if(inherits(study_boundary,"sf") && nrow(study_boundary)!=1){
    stop("study_boundary must only have a single record geometry")
  }
  if(inherits(study_boundary,"sfc") && length(study_boundary)!=1){
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

  # buff
  # check buff
  if(is.character(buffer)){
    buffer <- readr::parse_number(buffer)
  }
  if(dplyr::coalesce(buffer,0)>0){
    study_boundary <- sf::st_buffer(study_boundary, buffer)
  }

  # bbox
  if(bbox_of_aoi){
    study_boundary <- sf::st_bbox(study_boundary) %>%
      sf::st_as_sfc()
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
  if(
    !inherits(sf_data,"sf")
    && !inherits(sf_data,"sfc")
  ){stop("sf_data must be sf class object")}
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
# intermediate function 4
# function to make a Fortran topography file for QUIC-Fire O_o
#######################################################
quicfire_dtm_topofile <- function(
  dtm_rast
  , horizontal_resolution = 2
  , study_boundary = NULL
  , outdir = tempdir()
){
  # open dtm
  if(inherits(dtm_rast, "SpatRaster")) { # terra
    dtm_rast <- dtm_rast %>% terra::subset(1)
  }else if(inherits(dtm_rast,"RasterLayer")){ # raster
    dtm_rast <- terra::rast(dtm_rast) %>% terra::subset(1)
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
      dtm_rast <- terra::rast(fls[1]) %>% ## only reads one file...but what if multiple?
        terra::subset(1)

    }else if(any(stringr::str_ends(dtm_rast, ".*\\.(tif|tiff)$"))){
      # read it
      dtm_rast <- stringr::str_subset(dtm_rast, pattern = ".*\\.(tif|tiff)$")[1] %>%
        ## only reads one file...but what if multiple?
        terra::rast() %>%
        terra::subset(1)
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

  # check the study_boundary
  if(
    inherits(study_boundary,"sf")
    || inherits(study_boundary,"sfc")
  ){
    study_boundary <- get_custom_aoi(
      study_boundary = study_boundary
      , bbox_of_aoi = F, buffer = 0
      , reproject_epsg = terra::crs(dtm_rast) %>% sf::st_crs()
    )
  }

  # outdir
  outdir <- file.path(outdir)
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }

  # we need to get the raster in the desired resolution
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
  if(
    inherits(study_boundary,"sf")
    || inherits(study_boundary,"sfc")
  ){
    clipped_dtm <- terra::crop(dtm_rast, study_boundary)
  }
  #sum(is.na(values(clipped_dtm)))#count NANS

  #write the clipped dtm to a tif file just to be nice
  terra::writeRaster(
    x = clipped_dtm
    , filename = file.path(outdir,"dtm_Clipped.tif")
    , overwrite = TRUE
  ) #write clipped dtm raster to tif to check

  ### Write to FORTRAN File in format needed for QUIC-Fire (and FIRETEC)

  # Flip the raster over the y-axis
  raster_data <- terra::flip(clipped_dtm, direction = "vertical")

  # Extract the raster values
  values <- terra::values(raster_data) %>% c()

  # Open a connection to the unformatted Fortran file
  fortran_file <- file(file.path(outdir,"topo.dat"), "wb")

  # Write the data to the file
  writeBin(charToRaw("BRUH"), fortran_file) #Topo requires a header of 4 bits... this works and brings Sophie B. joy
  writeBin(values, fortran_file, size = 4)  # Assuming 4-byte (32-bit) floats

  # Close the file connection
  close(fortran_file)

  message(paste0(
    "exported QUIC-Fire topo.dat to:"
    , "\n ........ "
    , file.path(outdir,"topo.dat")
  ))

  return(clipped_dtm)
}


#######################################################
# intermediate function 5
# Format and export the data to work directly in TREES (default of 2m horizontal resolution) with Fuellist
#######################################################
export_to_TREES <- function(lines, data, box_coords, project_path,
                            litter=list( 'ilitter'  = 0,# 0 = no litter, 1 = litter
                                         'lrho'     = 4.667, #litter bulk density
                                         'lmoisture'= 0.06,  #litter moisture
                                         'lss'      = 0.0005,#litter sizescale
                                         'ldepth'   = 0.06), #litter depth
                            grass=list(  'igrass'   = 0,# 0 = no grass, 1 = grass
                                         'grho'     = 1.17, #grass bulk density
                                         'gmoisture'= 0.06,  #grass moisture
                                         'gss'      = 0.0005,#grass sizescale
                                         'gdepth'   = 0.27), #grass depth
                            topofile='flat',CBD_choice="cruz_tree_kg_per_m3",horizontal_resolution=2)
{


    data <- data %>% as.data.frame %>% dplyr::select(-geom)   #remove geometry column

    #Convert coordinates to be relative to SW corner of that geojson we made earlier
    data$x_coord = data$tree_x - box_coords$xmin
    data$y_coord = data$tree_y - box_coords$ymin

    #Rearrange the dataframe columns
    treelist <- data %>% dplyr::select(sp,x_coord,y_coord,tree_height_m,tree_cbh_m,crown_dia_m,max_crown_diam_height_m,CBD_choice,moist,ss) %>%
      mutate(across(where(is.numeric), ~ round(., 4))) #round to 4 decimal places for clarity (and TREES doesn't need that much precision)

    #remove trees with no data in row
    treelist <- na.omit(treelist)

    #write full George treelist to formatted txt file (MUST BE A TEXTFILE WITHOUT HEADERS)
    write.table(treelist, file=paste0(project_path,"Cloud2Trees_TreeList.txt"), row.names = FALSE,sep=" ",  col.names=FALSE)

    print("Exported Treelist for LANL TREES Program!")

    nz = ceiling(max(treelist$tree_height_m)) + 1 #max tree height + 1 m

    # Make Fuellist for TREES
    lines[4]  <- paste0("      nx  = "      ,box_coords$nx)#nx
    lines[5]  <- paste0("      ny  = "      ,box_coords$ny)#ny
    lines[6]  <- paste0("      nz  = "      ,nz)#nz
    lines[7]  <- paste0("      dx  = "      ,horizontal_resolution)
    lines[8]  <- paste0("      dy  = "      ,horizontal_resolution)
    lines[13] <- paste0("      topofile = " ,topofile) #This is always flat for QF, but can be applied as a pathfile to topo.dat for FIRETEC
    lines[36] <- paste0("      treefile = '",project_path,"Cloud2Trees_DomainSW.txt'")#treefile path
    lines[37] <- paste0("      ndatax = "   ,box_coords$width)#ndatax
    lines[38] <- paste0("      ndatay = "   ,box_coords$length)#ndatay

    #Litter
    lines[44] <- paste0("      ilitter = "  ,litter$ilitter)# ! Litter flag; 0=no litter, 1=basic litter, 2=DUET
    lines[47] <- paste0("      lrho = "     ,litter$lrho)
    lines[48] <- paste0("      lmoisture = ",litter$lmoisture)
    lines[49] <- paste0("      lss = "      ,litter$lss)
    lines[50] <- paste0("      ldepth = "   ,litter$ldepth)

    #Grass
    lines[71] <- paste0("      igrass = "   ,grass$igrass)
    lines[75] <- paste0("      grho = "     ,grass$grho)
    lines[76] <- paste0("      gmoisture = ",grass$gmoisture)
    lines[77] <- paste0("      gss = "      ,grass$gss)
    lines[78] <- paste0("      gdepth = "   ,grass$gdepth)

    print(lines)

    writeLines(lines, paste0(project_path,"fuellist"))

    print('Exported TREES Fuellist')

}
