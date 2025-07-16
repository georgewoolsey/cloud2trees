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
#' @param bbox_aoi logical. Should the study_boundary be transformed to a bounding box instead of it's original shape for determining the trees within the boundary?
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
  , bbox_aoi = F
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
      , bbox_aoi = dplyr::coalesce(bbox_aoi,FALSE)
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
      , aoi = study_boundary # will be the bbox as an sf if bbox_aoi
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
  , bbox_aoi = F
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
  if(bbox_aoi){
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
  fp <- file.path(normalizePath(outdir),"Lidar_Bounds.geojson")
  sf::st_write(
    bbox
    , fp
    , driver = "GeoJSON"
    , delete_dsn = TRUE
  )

  message(paste0(
    "exported QUIC-Fire domain to:"
    , "\n ........ "
    , fp
  ))

  # make box_coords df
  box_coords_df <- bbox %>%
    sf::st_as_sf() %>%
    sf::st_transform(sf::st_crs(sf_data)) %>%
    sf::st_set_geometry("geometry") %>%
    dplyr::mutate(
      xmin = xmin
      , ymin = ymin
      , nx = nx
      , ny= ny
      , width = width
      , length = length
    )

  return(list(
    quicfire_domain_df = box_coords_df
    , domain_path = fp
  ))

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
      search_dir_final_detected_ans <- search_dir_final_detected(file.path(dtm_rast))
      dtm_flist <- search_dir_final_detected_ans$dtm_flist

      # stop it if no files
      if(length(dtm_flist)<1){
        stop(paste0(
          "no DTM raster file found at: "
          , "\n .... "
          , file.path(normalizePath(dtm_rast))
        ))
      }
      # read it
      dtm_rast <- terra::rast(dtm_flist[1]) %>% ## only reads one file...but what if multiple?
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
      , bbox_aoi = F, buffer = 0
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
  fp <- file.path(normalizePath(outdir),"topo.dat")
  fortran_file <- file(fp, "wb")

  # Write the data to the file
  writeBin(charToRaw("BRUH"), fortran_file) #Topo requires a header of 4 bits... this works and brings Sophie B. joy
  writeBin(values, fortran_file, size = 4)  # Assuming 4-byte (32-bit) floats

  # Close the file connection
  close(fortran_file)

  message(paste0(
    "exported QUIC-Fire topo.dat to:"
    , "\n ........ "
    , fp
  ))

  return(list(
    dtm = clipped_dtm
    , topofile_path = fp
  ))
}

#######################################################
# intermediate function 5
# Read the fuellist example file into a vector of lines (use default version of file and then save to new location, keeping this version)
#######################################################
# Read the fuellist example file into a vector of lines (use default version of file and then save to new location, keeping this version)
# MAKE THIS FILE INHERENT IN THE PACKAGE... USE AS BASELINE FOR TREELIST
quicfire_get_fuellist <- function() {
  f <- system.file("extdata", "fuellist", package = "cloud2trees")
  lines <- readLines(f)
  return(lines)
}

#######################################################
# intermediate function 6
# check the fuellist
#######################################################
quicfire_check_fuellist <- function(fuellist) {
  # # get the parameters in order from fuellist
  #   quicfire_get_fuellist() %>%
  #     dplyr::as_tibble() %>%
  #     dplyr::rename(huh=1) %>%
  #     dplyr::mutate(
  #       huh = stringr::str_squish(huh)
  #       , row_number = dplyr::row_number()
  #     ) %>%
  #     # remove comments and blanks
  #     dplyr::filter(
  #       !stringr::str_starts(huh, pattern = "\\!")
  #       & huh!=""
  #       & !is.na(huh)
  #     ) %>%
  #     # get the name of the parameter
  #     dplyr::mutate(
  #       parameter = stringr::word(huh, sep = "=") %>% stringr::str_squish()
  #     ) %>%
  #     # dplyr::pull(parameter) %>%
  #     dplyr::pull(row_number) %>%
  #     paste(collapse = ",")


    # the parameters in order form fuellist
    orig <- dplyr::tibble(
      parameter = c(
        # firetech domain info
        "nx","ny","nz","dx","dy","dz","aa1","singlefuel","lreduced","topofile"
        # data import from existing files
        ,"ifuelin","inx","iny","inz","idx","idy","idz","iaa1","infuel","intopofile","rhoffile","ssfile","moistfile","afdfile"
        # input trees dataset info
        ,"itrees","tfuelbins","treefile","ndatax","ndatay","datalocx","datalocy"
        # litter switch
        ,"ilitter","litterconstant","lrho","lmoisture","lss","ldepth","windprofile","YearsSinceBurn","StepsPerYear"
        ,"relhum","grassstep","iFIA","FIA","randomwinds","litout","gmoistoverride","uavg","vavg","ustd","vstd"
        # grass switch
        ,"igrass","ngrass","grassconstant","grho","gmoisture","gss","gdepth"
        # option to output tree list and fuellist
        ,"verbose"
      )
      , row_number = c(
        # firetech domain info
        4,5,6,7,8,9,10,11,12,13
        # data import from existing files
        ,17,18,19,20,21,22,23,24,25,26,27,28,29,30
        # input trees dataset info
        ,34,35,36,37,38,39,40
        # litter switch
        ,44,46,47,48,49,50,53,54,55
        ,56,57,58,59,60,62,63,64,65,66,67
        # grass switch
        ,71,73,74,75,76,77,78
        # option to output tree list and fuellist
        ,82
      )
    )

    # check
    if(!inherits(fuellist,"character")){stop("this fuellist isn't even character")}
    # compare to orig
    comp <- orig %>%
      dplyr::anti_join(
        fuellist %>%
          dplyr::as_tibble() %>%
          dplyr::rename(huh=1) %>%
          dplyr::mutate(
            huh = stringr::str_squish(huh)
            , row_number = dplyr::row_number()
          ) %>%
          # remove comments and blanks
          dplyr::filter(
            !stringr::str_starts(huh, pattern = "\\!")
            & huh!=""
            & !is.na(huh)
          ) %>%
          # get the name of the parameter
          dplyr::mutate(
            parameter = stringr::word(huh, sep = "=") %>% stringr::str_squish()
          )
        , by = dplyr::join_by("parameter", "row_number")
      ) %>%
      dplyr::mutate(l = stringr::str_c(row_number,parameter,sep=":"))

    if(nrow(comp)!=0){
      stop(paste0(
        "fuellist has missing or unordered parameters listed below (expected_row:parameter). see `quicfire_get_fuellist()`"
        , "\n   "
        , paste(comp$l, collapse = ", ")
      ))
    }else{
      return(T)
    }
}
# quicfire_check_fuellist(quicfire_get_fuellist() %>% sample(66))

#######################################################
# intermediate function 87
# Format and export the data to work directly in TREES (default of 2m horizontal resolution) with Fuellist
# could eventually just use lanl trees to directly generate QF input?
# https://github.com/lanl/Trees/
#######################################################
make_lanl_trees_input <- function(
  tree_list
  , quicfire_domain_df
  , topofile = "flat" #This is always flat for QF, but can be applied as a pathfile to topo.dat (eg "C:\\data\\topo.dat") for FIRETEC
  , cbd_col_name = "landfire_tree_kg_per_m3" # "cruz_tree_kg_per_m3" "landfire_tree_kg_per_m3"
  , horizontal_resolution = 2
  , outdir = tempdir()
  # fuel settings
  , fuel_litter = list(
    "ilitter"  = 0 # 0 = no litter, 1 = litter
    , "lrho"     = 4.667 #litter bulk density
    , "lmoisture"= 0.06 #litter moisture
    , "lss"      = 0.0005 #litter sizescale
    , "ldepth"   = 0.06 #litter depth
  )
  , fuel_grass = list(
    "igrass"   = 0 # 0 = no grass, 1 = grass
    , "grho"     = 1.17 #grass bulk density
    , "gmoisture"= 0.06 #grass moisture
    , "gss"      = 0.0005 #grass sizescale
    , "gdepth"   = 0.27 #grass depth
  )
){
  # eventually have species-specific fuel settings
  # this is the parameter str
  fuel_trees <-
    dplyr::tibble(
      "sp" = c(1) #make species column where its all 1 (for now)...idk what the numeric codes in QF are, do you?
      , "moist" = c(1) #setting moisture at 100% here for now... can be a species-level input
      , "ss" = c(0.0005) #species-level input, putting as the common value for pine for now (doesn't actually matter for QF atm)
    ) %>%
      # this is all just QC that could be done on parameter
      dplyr::group_by(sp) %>%
      dplyr::filter(dplyr::row_number()==1) %>%
      dplyr::ungroup()

  ############!!!!!!!!!!!!!!! need to check str of fuel_litter and fuel_grass !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ############!!!!!!!!!!!!!!! need to check str of fuel_litter and fuel_grass !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ############!!!!!!!!!!!!!!! need to check str of fuel_litter and fuel_grass !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ############!!!!!!!!!!!!!!! need to check str of fuel_litter and fuel_grass !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ############!!!!!!!!!!!!!!! need to check str of fuel_litter and fuel_grass !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  # topofile
  if(
    stringr::str_ends(tolower(topofile),"\\.dat")
  ){
    # make it the proper path format
    topofile <- topofile %>%
      normalizePath() %>%
      stringr::str_replace_all(pattern = "/", replacement = "\\\\\\") %>% # hopefully TREES can read this filepath format?
      sQuote(q=F) # q=F > the undirectional ASCII quotation style is used
  }else if(tolower(topofile)!="flat"){
    stop("topofile must be either 'flat' or the path to a '.dat' file")
  }


  # outdir
  outdir <- file.path(outdir)
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }

  # Read the fuellist example file into a vector of lines (use default version of file and then save to new location, keeping this version)
  lines <- quicfire_get_fuellist()
  quicfire_check_fuellist(lines)

  # check the str of the quicfire_domain_df
  check_df_cols_all_missing(
    df = quicfire_domain_df
    , col_names = c(
      "xmin","ymin","nx","ny","width","length"
    )
    , check_vals_missing = T
  )

  # convert tree list to spatial points data
  tree_list <- check_spatial_points(tree_list)

  # check the str of the tree_list
  check_df_cols_all_missing(
    tree_list
    , col_names = c(
      "tree_x", "tree_y"
      , "tree_height_m", "tree_cbh_m", "crown_dia_m", "max_crown_diam_height_m"
      , cbd_col_name
    )
    , check_vals_missing = T
    , all_numeric = T
  )

  # set up data for return
  data <-
    tree_list %>%
    sf::st_drop_geometry() %>%
    #make species column where its all 1 (for now)...idk what the numeric codes in QF are, do you?
    dplyr::mutate(sp = 1) %>% ### eventually this can be the QF species codes
    # add fuel_trees
    dplyr::left_join(fuel_trees, by = "sp") %>%
    #Convert coordinates to be relative to SW corner of that geojson we made earlier
    dplyr::mutate(
      x_coord = tree_x - quicfire_domain_df$xmin[1]
      , y_coord = tree_y - quicfire_domain_df$ymin[1]
    ) %>%
    # select and rearrange the columns
    dplyr::select(dplyr::all_of(c(
      "sp", "x_coord", "y_coord"
      , "tree_height_m", "tree_cbh_m", "crown_dia_m", "max_crown_diam_height_m"
      , cbd_col_name
      , "moist", "ss"
    ))) %>%
    #round to 4 decimal places for clarity (and TREES doesn't need that much precision)
    dplyr::mutate(dplyr::across(
      dplyr::where(is.numeric)
      , ~ round(.x, 4))
    ) %>%
    #remove trees with no data in row
    tidyr::drop_na()  # remove rows with any NA values in any column
    # # tidyr::drop_na(col1,col2) # remove rows with NA values in select columns
    # # dplyr::filter( # remove rows with NA values in all columns...i.e. keep rows where AT LEAST ONE column is NOT NA
    # #   dplyr::if_any(
    # #     dplyr::everything()
    # #     , ~ !is.na(.x))
    # # )

  # stop it if too many trees dropped
  if( nrow(data) < floor(nrow(tree_list)*0.9) ){
    stop("more than 10% of data dropped due to missing values in tree attributes...fill in missing data first")
  }

  #write full cloud2trees treelist to formatted txt file (MUST BE A TEXTFILE WITHOUT HEADERS)
  treefile_path <- file.path(normalizePath(outdir),"Cloud2Trees_TreeList.txt") %>%
    stringr::str_replace_all(pattern = "/", replacement = "\\\\\\") # hopefully TREES can read this filepath format?
  write.table(
    data
    , file = treefile_path
    , row.names = FALSE
    , sep = " "
    , col.names = FALSE
  )

  message(paste0(
    "exported Treelist for LANL TREES program to:"
    , "\n ........ "
    , treefile_path
  ))

  # enclose the fpath string in single quotes
  treefile_path_q <- sQuote(treefile_path, q = F) # q=F > the undirectional ASCII quotation style is used

  ###############################
  #### !!!! now the fuellist
  ## could make this a separate fn eventually
  ###############################
  nz <- ceiling(max(data$tree_height_m)) + 1 #max tree height + 1 m

  # Make Fuellist for TREES
  lines[4]  <- paste0("      nx  = "      ,quicfire_domain_df$nx[1])#nx
  lines[5]  <- paste0("      ny  = "      ,quicfire_domain_df$ny[1])#ny
  lines[6]  <- paste0("      nz  = "      ,nz)#nz
  lines[7]  <- paste0("      dx  = "      ,horizontal_resolution)
  lines[8]  <- paste0("      dy  = "      ,horizontal_resolution)
  lines[13] <- paste0("      topofile = " ,topofile) #This is always flat for QF, but can be applied as a pathfile to topo.dat for FIRETEC
  lines[36] <- paste0("      treefile = " ,treefile_path_q)#treefile path
  lines[37] <- paste0("      ndatax = "   ,quicfire_domain_df$width[1])#ndatax
  lines[38] <- paste0("      ndatay = "   ,quicfire_domain_df$length[1])#ndatay

  #Litter
  lines[44] <- paste0("      ilitter = "  ,fuel_litter$ilitter[1])# ! Litter flag; 0=no litter, 1=basic litter, 2=DUET
  lines[47] <- paste0("      lrho = "     ,fuel_litter$lrho[1])
  lines[48] <- paste0("      lmoisture = ",fuel_litter$lmoisture[1])
  lines[49] <- paste0("      lss = "      ,fuel_litter$lss[1])
  lines[50] <- paste0("      ldepth = "   ,fuel_litter$ldepth[1])

  #Grass
  lines[71] <- paste0("      igrass = "   ,fuel_grass$igrass[1])
  lines[75] <- paste0("      grho = "     ,fuel_grass$grho[1])
  lines[76] <- paste0("      gmoisture = ",fuel_grass$gmoisture[1])
  lines[77] <- paste0("      gss = "      ,fuel_grass$gss[1])
  lines[78] <- paste0("      gdepth = "   ,fuel_grass$gdepth[1])

  # print(lines)
  fuellist_path <- file.path(normalizePath(outdir),"fuellist")
  writeLines(lines, fuellist_path)

  message(paste0(
    "exported TREES program fuellist to:"
    , "\n ........ "
    , fuellist_path
  ))

  return(list(
    treelist = data
    , fuellist_path = fuellist_path
    , treelist_path = treefile_path
  ))
}






# library(tidyverse)
# # customize the aoi settings and clip the tree list
# clip_tree_list_aoi_ans <- clip_tree_list_aoi(
#   tree_list = sf::st_read("c:/data/usfs/dod_cloud2trees_demo/data/SycanMarsh/als_2021_processing/point_cloud_processing_delivery/final_detected_tree_tops.gpkg")
#   , study_boundary = sf::st_read("c:/data/usfs/dod_cloud2trees_demo/data/SycanMarsh/Sycan_2A.shp")
#   , bbox_aoi = T
#   , buffer = 50
#   , reproject_epsg = NULL
# )
# clip_tree_list_aoi_ans %>% names()
# inherits(clip_tree_list_aoi_ans$aoi,"sfc")
# sf::st_crs(clip_tree_list_aoi_ans$aoi)
# ggplot() + geom_sf(data=clip_tree_list_aoi_ans$aoi)
# # do the domain thing
# quicfire_domain_df <- quicfire_define_domain(
#   sf_data = clip_tree_list_aoi_ans$aoi
#   , horizontal_resolution = 2
#   , outdir = tempdir()
# )
# quicfire_domain_df %>% dplyr::glimpse()
# quicfire_domain_df %>% sf::st_crs(parameters = T)
# nrow(quicfire_domain_df)
# # don't forget the fortran flipped over rast
# quicfire_dtm_topofile_ans <- quicfire_dtm_topofile(
#   dtm_rast = terra::rast("c:/data/usfs/dod_cloud2trees_demo/data/SycanMarsh/als_2021_processing/point_cloud_processing_delivery/dtm_1m.tif")
#   , horizontal_resolution = 2
#   , study_boundary = quicfire_domain_df
#   , outdir = tempdir()
# )
# # terra::plot(quicfire_dtm_topofile_ans)
# # terra::plot(
# #   quicfire_domain_df %>%
# #     sf::st_transform(terra::crs(quicfire_dtm_topofile_ans)) %>%
# #     terra::vect()
# #   , add = T, border = "blue", col = NA
# #   , lwd = 22
# # )
#
# quicfire_get_fuellist() %>% length()
