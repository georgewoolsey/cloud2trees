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
#' If no boundary given, the AOI will be built from location of trees in the tree list.
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
  if(!inherits(study_boundary,"sf")){stop("study_boundary must be sf class object")}
  if(is.na(sf::st_crs(tree_list))){stop("tree_list does not have a CRS")}
  if(is.na(sf::st_crs(study_boundary))){stop("study_boundary does not have a CRS")}
  
  # check epsg buff
  if(is.character(reproject_epsg)){
    reproject_epsg <- readr::parse_number(reproject_epsg) 
  }
  if(is.character(buffer)){
    buffer <- readr::parse_number(buffer) 
  }
  
  # bounds
  if(
    !all( sf::st_is(study_boundary, c("POLYGON","MULTIPOLYGON")) )
  ){
    stop("study_boundary must contain POLYGON type geometry only")
  }
  if(nrow(study_boundary)!=1){
    stop("study_boundary must only have a single record geometry")
  }

  ###############################
  # data transformations
  ###############################
  # convert tree list to spatial points data
  tree_tops <- check_spatial_points(tree_list, crs)
  
  # reproj
  if(!is.null(reproject_epsg) && is.numeric(reproject_epsg)){
    tree_tops <- tree_tops %>% sf::st_transform(crs = reproject_epsg)
  }
  
  # bbox
  if(bbox_of_aoi){
    study_boundary <- sf::st_bbox(study_boundary) %>% 
      sf::st_as_sfc() %>% 
      sf::st_transform(sf::st_crs(tree_tops))
  }else{
    study_boundary <- study_boundary %>% sf::st_transform(sf::st_crs(tree_tops))
  }
  
  # buff
  if(dplyr::coalesce(buffer,0)>0){
    study_boundary <- sf::st_buffer(study_boundary, buffer)
  }
  
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
    warning("no tree_list found within study_boundary bounds")
    return(ret)
  }else{
    return(ret)
  }

}
#######################################################
# intermediate function 2
#######################################################
