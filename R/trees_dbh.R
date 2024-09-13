#' @title Estimate DBH for a tree list based on height
#'
#' @description
#' `trees_dbh()` uses the input tree list (e.g. as exported by [raster2trees()]) with the columns
#' `treeID`, `tree_x`, `tree_y`, and `tree_height_m` to estimate tree DBH.
#'
#' The process using the local model (default):
#'
#' * Use the TreeMap (see reference) FIA plot data in the area of the tree list to estimate the height-DBH allometry relationship
#' * Use the model estimated above to predict DBH based on tree height
#'
#' @param tree_list data.frame. A data frame with the columns `treeID`, `tree_x`, `tree_y`, and `tree_height_m`.
#' If an `sf` class object with POINT geometry (see [sf::st_geometry_type()]), the program will use the data "as-is" and only require the `treeID` and `tree_height_m` columns.
#' @param crs string. A crs string as returned from [sf::st_crs()] or the EPSG code of the x,y coordinates.
#' Defaults to the crs of the `tree_list` data if of class "sf".
#' @param study_boundary sf. If you want to scale per ha calculations, provide the geography of the study boundary
#'
#' @references
#' [https://doi.org/10.2737/RDS-2021-0074](https://doi.org/10.2737/RDS-2021-0074)
#' Riley, Karin L.; Grenfell, Isaac C.; Finney, Mark A.; Shaw, John D. 2021. TreeMap 2016: A tree-level model of the forests of the conterminous United States circa 2016. Fort Collins, CO: Forest Service Research Data Archive.
#'
#' @return Returns a spatial data frame of individual trees.
#'
#' @examples
#'  \dontrun{
#'  o <- "../data"
#'  i <- "../data/lasdata"
#'  r <- cloud2trees::cloud2raster(output_dir = o, input_las_dir = i)
#'  r %>% names()
#'  r$dtm_rast %>% terra::plot()
#'  r$chm_rast %>% terra::plot()
#'  r$create_project_structure_ans %>% dplyr::glimpse()
#'  r$chunk_las_catalog_ans$process_data %>% dplyr::glimpse()
#'  r$chunk_las_catalog_ans$grid_subset_switch
#'  r$chunk_las_catalog_ans$las_ctg@data %>% dplyr::glimpse()
#'  r$normalize_flist
#'  }
#' @export
#'
trees_dbh <- function(
  tree_list
  , crs = NA
  , study_boundary = NA
) {
  ##################################
  # convert to spatial points data
  ##################################
  if(inherits(tree_list, "sf")){
    # if points, just use it
    if( min(sf::st_is(tree_list, type = c("POINT", "MULTIPOINT"))) == 1 ){
      tree_tops <- tree_list
    }else{ # if spatial but not points, drop geom and set to points
      tree_tops <- tree_list %>%
        sf::st_drop_geometry() %>%
        sf::st_as_sf(
          coords = c("tree_x", "tree_y"), crs = sf::st_crs(tree_list)
          , remove = F
        )
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
      )
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
        , max_tree_height_m = max(comp_tree_height_m)
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
