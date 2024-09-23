#' @title Calculate competition metrics for a tree list
#'
#' @description
#' `trees_competition()` uses the input tree list (e.g. as exported by [raster2trees()]) with the columns
#' `treeID`, `tree_x`, `tree_y`, and `tree_height_m` to calculate competition metrics at the tree level.
#'
#' Competition metrics returned include:
#'
#' * Distance to the nearest neighbor (`comp_dist_to_nearest_m`)
#' * Trees per ha within a 5m radius (`comp_trees_per_ha`)
#' * The relative tree height (`comp_relative_tree_height` = height_tree / height_max) within a 5m radius where a value of 1 indicates the tallest tree.
#'
#' @param tree_list data.frame. A data frame with the columns `treeID`, `tree_x`, `tree_y`, and `tree_height_m`.
#' If an `sf` class object with POINT geometry (see [sf::st_geometry_type()]), the program will use the data "as-is" and only require the `treeID` and `tree_height_m` columns.
#' @param crs string. A crs string as returned from [sf::st_crs()] or the EPSG code of the x,y coordinates.
#' Defaults to the crs of the `tree_list` data if of class "sf".
#' @param competition_buffer_m number. Set buffer around tree (m) to calculate competition metrics
#' @param study_boundary sf. If you want to scale per ha calculations, provide the geography of the study boundary
#' @param search_dist_max number. Maximum search distance (m) to nearest tree. Larger search distances will increase processing time and possibly result in memory issues.
#' If no competition trees are found within this distance, the return column `comp_dist_to_nearest_m` = `search_dist_max` parameter.
#'
#' @references
#' [https://doi.org/10.3390/f13122077](https://doi.org/10.3390/f13122077)
#' Tinkham et al. (2022). Modeling the missing DBHs: Influence of model form on UAV DBH characterization. Forests, 13(12), 2077.
#'
#' @return Returns a spatial data frame of individual trees.
#'
#' @examples
#'  \dontrun{
#'  # example tree list
#'  tl <- dplyr::tibble(
#'      treeID = c(1:21)
#'      , tree_x = rnorm(n=21, mean = 458064, sd = 11)
#'      , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
#'      , tree_height_m = exp(rgamma(n = 21, shape = (7/4)^2, rate = (4^2)/7))
#'    )
#'  # call the function
#'  tl_comp <- trees_competition(tree_list = tl, crs = "32613")
#'  # what?
#'  tl_comp %>% class()
#'  tl_comp %>% dplyr::select(tidyselect::starts_with("comp_")) %>% dplyr::glimpse()
#'  tl_comp %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=comp_dist_to_nearest_m))
#'  }
#' @export
#'
trees_competition <- function(
  tree_list
  , crs = NA
  , competition_buffer_m = 5
  , study_boundary = NA
  , search_dist_max = 10
) {
  ##################################
  # ensure that tree height data exists
  ##################################
  f <- tree_list %>% names()
  if(length(f)==0){f <- ""}
  if(
    max(grepl("tree_height_m", f))==0
  ){
    stop(paste0(
      "`tree_list` data must contain `tree_height_m` column to estimate DBH."
      , "\nRename the height column if it exists and ensure it is in meters."
    ))
  }
  if(
    max(grepl("treeID", f))==0
  ){
    stop(paste0(
      "`tree_list` data must contain `treeID` column to estimate DBH."
      , "\nProvide the `treeID` as a unique identifier of individual trees."
    ))
  }

  ##################################
  # convert to spatial points data
  ##################################
  if(inherits(tree_list, "sf")){
    # if points, just use it
    if( min(sf::st_is(tree_list, type = c("POINT", "MULTIPOINT"))) == 1 ){
      tree_tops <- tree_list %>%
        dplyr::mutate(treeID = as.character(treeID), tree_height_m = as.numeric(tree_height_m))
    }else{ # if spatial but not points, drop geom and set to points
      tree_tops <- tree_list %>%
        sf::st_drop_geometry() %>%
        sf::st_as_sf(
          coords = c("tree_x", "tree_y"), crs = sf::st_crs(tree_list)
          , remove = F
        ) %>%
        dplyr::mutate(treeID = as.character(treeID), tree_height_m = as.numeric(tree_height_m))
    }
  }else{ # not spatial data
    # convert from data.frame to spatial points
    if(!inherits(tree_list, "data.frame")){
      stop("must pass a data.frame or sf object to the `tree_list` parameter")
    }
    if(is.na(crs) | is.na(readr::parse_number(as.character(crs)))){
      stop("must provide the EPSG code in `crs` parameter for the projection of x,y data")
    }
    tree_tops <- tree_list %>%
      sf::st_as_sf(
        coords = c("tree_x", "tree_y")
        , crs = paste0( "EPSG:", readr::parse_number(as.character(crs)) )
        , remove = F
      ) %>%
      dplyr::mutate(treeID = as.character(treeID), tree_height_m = as.numeric(tree_height_m))
  }

  # check for duplicate treeID
  if(
    nrow(tree_tops) != length(unique(tree_tops$treeID))
  ){
    stop("Duplicates found in the treeID column. Please remove duplicates and try again.")
  }

  ####################################################################
  # Calculate local tree competition metrics
  ####################################################################
  if(inherits(study_boundary, "sf") | inherits(study_boundary, "sfc")){
    ### how much of the buffered tree area is within the study boundary?
      # use this to scale the TPA estimates below
    tree_tops_pct_buffer_temp <- tree_tops %>%
      # buffer point
      sf::st_buffer(competition_buffer_m) %>%
      dplyr::mutate(
        point_buffer_area_m2 = as.numeric(sf::st_area(.))
      ) %>%
      # intersect with study bounds
      sf::st_intersection(
        study_boundary %>%
          sf::st_union() %>%
          sf::st_as_sf() %>%
          sf::st_transform(sf::st_crs(tree_tops))
      ) %>%
      # calculate area of buffer within study
      dplyr::mutate(
        buffer_area_in_study_m2 = as.numeric(sf::st_area(.))
      ) %>%
      sf::st_drop_geometry() %>%
      dplyr::select(treeID, buffer_area_in_study_m2)
  }else{
    tree_tops_pct_buffer_temp <- dplyr::tibble(treeID = character(0), buffer_area_in_study_m2 = numeric(0))
  }
    ### use the tree top location points to get competition metrics
    comp_tree_tops_temp <- tree_tops %>%
      # buffer point
      sf::st_buffer(competition_buffer_m) %>%
      dplyr::select(treeID, tree_height_m) %>%
      # spatial join with all tree points
      sf::st_join(
        tree_tops %>%
          dplyr::select(treeID, tree_height_m) %>%
          dplyr::rename_with(
            .fn = ~ paste0("comp_",.x)
            , .cols = tidyselect::everything()[
                -dplyr::any_of(c("geometry"))
              ]
          )
      ) %>%
      sf::st_drop_geometry() %>%
      # calculate metrics by treeID
      dplyr::group_by(treeID,tree_height_m) %>%
      dplyr::summarise(
        n_trees = dplyr::n()
        , max_tree_height_m = max(comp_tree_height_m, na.rm = T)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(
        tree_tops_pct_buffer_temp
        , by = dplyr::join_by("treeID")
      ) %>%
      dplyr::mutate(
        comp_trees_per_ha = (n_trees/dplyr::coalesce(buffer_area_in_study_m2,1))*10000
        , comp_relative_tree_height = tree_height_m/max_tree_height_m
      ) %>%
      dplyr::select(
        treeID, comp_trees_per_ha, comp_relative_tree_height
      )

    ### calculate distance to nearest neighbor
      # get trees within radius
      dist_tree_tops_temp <- tree_tops %>%
        dplyr::select(treeID) %>%
        # buffer point
        sf::st_buffer(search_dist_max) %>%
        # spatial join with all tree points
        sf::st_join(
          tree_tops %>%
            dplyr::select(treeID, tree_x, tree_y) %>%
            dplyr::rename(treeID2=treeID)
        ) %>%
        dplyr::filter(treeID != treeID2)

      # calculate row by row distances
      dist_tree_tops_temp <- dist_tree_tops_temp %>%
        sf::st_centroid() %>%
        sf::st_set_geometry("geom1") %>%
        dplyr::bind_cols(
          dist_tree_tops_temp %>%
            sf::st_drop_geometry() %>%
            dplyr::select("tree_x", "tree_y") %>%
            sf::st_as_sf(coords = c("tree_x", "tree_y"), crs = sf::st_crs(tree_tops)) %>%
            sf::st_set_geometry("geom2")
        ) %>%
        dplyr::mutate(
          comp_dist_to_nearest_m = sf::st_distance(geom1, geom2, by_element = T) %>% as.numeric()
        ) %>%
        sf::st_drop_geometry() %>%
        dplyr::select(treeID,comp_dist_to_nearest_m) %>%
        dplyr::group_by(treeID) %>%
        dplyr::summarise(comp_dist_to_nearest_m = min(comp_dist_to_nearest_m, na.rm = T)) %>%
        dplyr::ungroup()

    ### join with original tree tops data
    tree_tops <- tree_tops %>%
      dplyr::left_join(
        comp_tree_tops_temp
        , by = dplyr::join_by("treeID")
      ) %>%
      # add distance
      dplyr::left_join(
        dist_tree_tops_temp
        , by = dplyr::join_by("treeID")
      ) %>%
      dplyr::mutate(comp_dist_to_nearest_m = dplyr::coalesce(comp_dist_to_nearest_m,search_dist_max))

  # return
  return(tree_tops)
}
