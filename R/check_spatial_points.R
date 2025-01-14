#' @title Check a data frame for spatial point data. Convert to points if needed.
#'
#' @description
#' Check a data frame for spatial point data. Convert to points if needed.
#'
#' @param tree_list data.frame. A data frame with the columns `treeID`, `tree_x`, `tree_y`, and `tree_height_m`.
#' If an `sf` class object with POINT geometry (see [sf::st_geometry_type()]), the program will use the data "as-is" and only require the `treeID` and `tree_height_m` columns.
#' @param crs string. A crs string as returned from [sf::st_crs()] or the EPSG code of the x,y coordinates.
#' Defaults to the crs of the `tree_list` data if of class "sf".
#'
#' @keywords internal
#'
check_spatial_points <- function(
  tree_list
  , crs = NA
){
  ##################################
  # ensure that tree data exists
  ##################################
  nms <- tree_list %>% names() %>% dplyr::coalesce("")
  if(
    !(stringr::str_equal(nms, "treeID") %>% any())
  ){
    stop(paste0(
      "`tree_list` data must contain `treeID` column."
      , "\nProvide the `treeID` as a unique identifier of individual trees."
    ))
  }

  # check for duplicate treeID
  if(
    nrow(tree_list) != length(unique(tree_list$treeID))
  ){
    stop("Duplicates found in the treeID column. Please remove duplicates and try again.")
  }

  ##################################
  # convert to spatial points data
  ##################################
  if(inherits(tree_list, "sf")){
    # if points, just use it
    if( sf::st_is(tree_list, type = c("POINT", "MULTIPOINT")) %>% all() ){
      tree_tops <- tree_list %>%
        dplyr::mutate(
          tree_x = sf::st_coordinates(.)[,1]
          , tree_y = sf::st_coordinates(.)[,2]
        ) %>%
        dplyr::mutate(treeID = as.character(treeID))
    }else{ # if spatial but not points, drop geom and set to points
      if(
        !(stringr::str_equal(nms, "tree_x") %>% any())
        && !(stringr::str_equal(nms, "tree_y") %>% any())
      ){ # doesn't contain x,y
        tree_tops <- tree_list %>%
          sf::st_centroid() %>%
          dplyr::mutate(
            tree_x = sf::st_coordinates(.)[,1]
            , tree_y = sf::st_coordinates(.)[,2]
          )
      }
      tree_tops <- tree_list %>%
        sf::st_drop_geometry() %>%
        sf::st_as_sf(
          coords = c("tree_x", "tree_y"), crs = sf::st_crs(tree_list)
          , remove = F
        ) %>%
        dplyr::mutate(treeID = as.character(treeID))
    }
  }else{ # not spatial data
    # convert from data.frame to spatial points
    if(!inherits(tree_list, "data.frame")){
      stop("must pass a data.frame or sf object to the `tree_list` parameter")
    }
    if(is.na(crs) || is.na(readr::parse_number(as.character(crs)))){
      stop("must provide the EPSG code in `crs` parameter for the projection of x,y data")
    }
    if(
      !(stringr::str_equal(nms, "tree_x") %>% any())
      || !(stringr::str_equal(nms, "tree_y") %>% any())
    ){ # doesn't contain x,y
      stop("must provide the columns `tree_x` and `tree_y`")
    }
    tree_tops <- tree_list %>%
      sf::st_as_sf(
        coords = c("tree_x", "tree_y")
        , crs = paste0( "EPSG:", readr::parse_number(as.character(crs)) )
        , remove = F
      ) %>%
      dplyr::mutate(treeID = as.character(treeID))
  }

  return(tree_tops)
}
