#' @title Use raw .las|.laz files to generate CHM, DTM, and a tree list
#'
#' @description
#' `cloud2trees()` is an all-in-one function to process raw .las|.laz files
#' to generate a CHM raster (.tif), a DTM raster (.tif), and a tree list with tree location, height, and DBH
#' The order of operations is:
#'
#' * Generate a CHM from the point cloud using [cloud2raster()]
#' * Perform individual tree detection using [raster2trees()]
#' * Quantify individual tree competition metrics using [trees_competition()] (*if set to TRUE*)
#' * Extract tree DBH values from the normalized point cloud using [treels_stem_dbh()] (*if set to TRUE*)
#' * Model tree DBH values using [trees_dbh()] (*if set to TRUE*)
#'
#' See the documentation for each individual function called for more details.
#'
#' @inheritParams cloud2raster
#' @param ws numeric or function. Length or diameter of the moving window used to detect the local
#' maxima in the units of the input data (usually meters). If it is numeric a fixed window size is used.
#' If it is a function, the function determines the size of the window at any given location on the canopy.
#' By default function takes the height of a given pixel as its only argument and return the
#' desired size of the search window when centered on that pixel.
#' @param estimate_tree_dbh logical. Should tree DBH be estimated? See [trees_dbh()].
#' @param max_dbh numeric. Set the largest tree diameter (m) expected in the point cloud
#' @param dbh_model string. Set the model to use for local dbh-height allometry. Can be "rf" for random forest or "lin" for linear
#' @param estimate_dbh_from_cloud logical. Should DBH be estimated from the point cloud? See [treels_stem_dbh()]. Setting to `TRUE` may significantly increase processing time.
#' @param estimate_tree_competition logical. Should tree competition metrics be calculated? See [trees_competition()]. Setting to `TRUE` may slightly increase processing time.
#' @param competition_buffer_m number. Set buffer around tree (m) to calculate competition metrics
#' @param search_dist_max number. Maximum search distance (m) to nearest tree for competition. Larger search distances will increase processing time and possibly result in memory issues.
#' If no competition trees are found within this distance, the return column `comp_dist_to_nearest_m` = `search_dist_max` parameter.
#'
#'
#' @return Returns the goods.
#' Exports files of the goods to new folders "point_cloud_processing_delivery" and "point_cloud_processing_temp" in the
#' `output_dir` defined by the user in the function call.
#'
#' @examples
#'  \dontrun{
#'  o <- "hey"
#'  }
#' @export
#'
cloud2trees <- function(
  output_dir
  , input_las_dir
  , input_treemap_dir = paste0(system.file(package = "cloud2trees"),"/extdata/treemap")
  , accuracy_level = 2
  , max_ctg_pts = 70e6
  , max_area_m2 = 90e6
  , transform = FALSE
  , new_crs = NA
  , old_crs = NA
  , keep_intrmdt = FALSE
  , dtm_res_m = 1
  , chm_res_m = 0.25
  , min_height = 2
  , max_height = 70
  , ws = function(x){
         y <- dplyr::case_when(
           is.na(x) ~ 0.001
           , x < 0 ~ 0.001
           , x < 2 ~ 1
           , x > 30 ~ 5
           , TRUE ~ 0.75 + (x * 0.14)
          )
         return(y)
    }
  , estimate_tree_dbh = FALSE
  , max_dbh = 2
  , dbh_model = "lin"
  , estimate_dbh_from_cloud = FALSE
  , estimate_tree_competition = FALSE
  , competition_buffer_m = 5
  , search_dist_max = 10
  , overwrite = TRUE
){
  ####################################################################
  # cloud2trees::cloud2raster()
  ####################################################################
  cloud2raster_ans <- cloud2trees::cloud2raster(
      output_dir = output_dir
      , input_las_dir = input_las_dir
      , input_treemap_dir = input_treemap_dir
      , accuracy_level = accuracy_level
      , max_ctg_pts = max_ctg_pts
      , max_area_m2 = max_area_m2
      , transform = transform
      , new_crs = new_crs
      , old_crs = old_crs
      , keep_intrmdt = keep_intrmdt
      , dtm_res_m = dtm_res_m
      , chm_res_m = chm_res_m
      , min_height = min_height
      , max_height = max_height
      , overwrite = overwrite
  )
  # cloud2raster_ans %>% names()
  # cloud2raster_ans$dtm_rast %>% terra::plot()
  # cloud2raster_ans$chm_rast %>% terra::plot()
  # cloud2raster_ans$create_project_structure_ans %>% dplyr::glimpse()
  # cloud2raster_ans$chunk_las_catalog_ans$process_data %>% dplyr::glimpse()
  # cloud2raster_ans$chunk_las_catalog_ans$is_chunked_grid
  # cloud2raster_ans$chunk_las_catalog_ans$las_ctg@data %>% dplyr::glimpse()
  # cloud2raster_ans$normalize_flist

  ####################################################################
  # cloud2trees::raster2trees()
  ####################################################################
  raster2trees_ans <- raster2trees(
    chm_rast = cloud2raster_ans$chm_rast
    , outfolder = cloud2raster_ans$create_project_structure_ans$delivery_dir
    , ws = ws
    , min_height = min_height
    , tempdir = cloud2raster_ans$create_project_structure_ans$temp_dir
  )

  # raster2trees_ans %>% class()
  # raster2trees_ans %>% typeof()
  # raster2trees_ans %>% dplyr::glimpse()
  # raster2trees_ans %>% sf::st_geometry_type() %>% table()
  # raster2trees_ans %>% sf::st_is(type = c("POLYGON", "MULTIPOLYGON")) %>% table()
  # raster2trees_ans %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(fill=tree_height_m))
  #
  # raster2trees_ans %>% ggplot2::ggplot() +
  #   ggplot2::geom_tile(
  #     data = cloud2raster_ans$chm_rast %>%
  #       as.data.frame(xy = T) %>%
  #       dplyr::rename(f=3)
  #     , mapping = ggplot2::aes(x=x,y=y,fill=f)
  #     , color = NA
  #   ) +
  #   ggplot2::geom_sf(ggplot2::aes(color=tree_height_m), fill = NA, lwd = 1.2) +
  #   ggplot2::scale_fill_viridis_c(option = "plasma") +
  #   ggplot2::scale_color_distiller(palette = "Greys", direction = 1) +
  #   ggplot2::theme_void() +
  #   ggplot2::theme(legend.position = "none")

  ####################################################################
  # cloud2trees::treels_stem_dbh()
  ####################################################################
  if(estimate_dbh_from_cloud==T){
    # treels stuff
    treels_dbh_locations <- cloud2trees::treels_stem_dbh(
      folder = cloud2raster_ans$normalize_flist
      , outfolder = cloud2raster_ans$create_project_structure_ans$treels_dir
      , min_height = min_height
      , max_dbh = max_dbh
      , chunk_these = !cloud2raster_ans$chunk_las_catalog_ans$is_chunked_grid
    )

    # treels_dbh_locations %>% class()
    # sf::st_geometry_type(treels_dbh_locations)
    # treels_dbh_locations %>% sf::st_is(type = c("POINT", "MULTIPOINT")) %>% min()
    # treels_dbh_locations %>% typeof()
    # treels_dbh_locations %>% dplyr::glimpse()
    # treels_dbh_locations %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(size = dbh_cm))
  }

  ####################################################################
  # cloud2trees::trees_competition()
  ####################################################################
  if(estimate_tree_competition==T){
    # trees_competition
    trees_competition_ans <- trees_competition(
      tree_list = raster2trees_ans
      , competition_buffer_m = competition_buffer_m
      , study_boundary = cloud2raster_ans$chunk_las_catalog_ans$las_ctg@data$geometry
      , search_dist_max = search_dist_max
    )

    # trees_competition_ans %>% class()
    # trees_competition_ans %>% dplyr::glimpse()
    # trees_competition_ans %>% dplyr::select(tidyselect::starts_with("comp_")) %>% dplyr::glimpse()
    # trees_competition_ans %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=tree_height_m))
    # trees_competition_ans %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=comp_dist_to_nearest_m))
  }

  ####################################################################
  # cloud2trees::trees_dbh()
  ####################################################################
  if(estimate_tree_dbh==T & estimate_dbh_from_cloud==T){
    # trees_dbh
    trees_dbh_ans <- trees_dbh(
      tree_list = raster2trees_ans
      , study_boundary = cloud2raster_ans$chunk_las_catalog_ans$las_ctg@data$geometry
      , dbh_model = dbh_model
      , treels_dbh_locations = treels_dbh_locations
      , input_treemap_dir = cloud2raster_ans$create_project_structure_ans$input_treemap_dir
      , outfolder = cloud2raster_ans$create_project_structure_ans$delivery_dir
    )
  }else if(estimate_tree_dbh==T){
    # trees_dbh
    trees_dbh_ans <- trees_dbh(
      tree_list = raster2trees_ans
      , study_boundary = cloud2raster_ans$chunk_las_catalog_ans$las_ctg@data$geometry
      , dbh_model = dbh_model
      , treels_dbh_locations = NA
      , input_treemap_dir = cloud2raster_ans$create_project_structure_ans$input_treemap_dir
      , outfolder = cloud2raster_ans$create_project_structure_ans$delivery_dir
    )
  }

  # trees_dbh_ans %>% class()
  # trees_dbh_ans %>% dplyr::glimpse()
  # trees_dbh_ans %>% dplyr::select(tidyselect::contains("dbh_cm")) %>% dplyr::glimpse()
  # trees_dbh_ans %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=dbh_cm))
  # trees_dbh_ans %>% ggplot2::ggplot(ggplot2::aes(x = tree_height_m, y = dbh_cm)) + ggplot2::geom_point()

  ####################################################################
  # return
  ####################################################################
    # return
    return(list(
      cloud2raster_ans = cloud2raster_ans
      , raster2trees_ans = raster2trees_ans
      # dtm_rast = dtm_rast
      # , chm_rast = chm_rast
      # , create_project_structure_ans = config
      # , chunk_las_catalog_ans = chunk_las_catalog_ans
      # , normalize_flist = normalize_flist
    ))
}
