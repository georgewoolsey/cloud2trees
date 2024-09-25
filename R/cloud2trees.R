#' @title Use raw .las|.laz files to generate CHM, DTM, and a tree list
#'
#' @description
#' `cloud2trees()` is an all-in-one function to process raw .las|.laz files
#' to generate a CHM raster (.tif), a DTM raster (.tif), and a tree list with tree location, height, and DBH.
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
#'  # test las file but this could also be a directory path with >1 .las|.laz files
#'  i <- system.file("extdata", "MixedConifer.laz", package="lidR")
#'  # run it
#'  cloud2trees_ans <- cloud2trees::cloud2trees(output_dir = tempdir(), input_las_dir = i)
#'  # what is it?
#'  cloud2trees_ans %>% names()
#'  # there's a DTM
#'  cloud2trees_ans$dtm_rast %>% terra::plot()
#'  # there's a CHM
#'  cloud2trees_ans$chm_rast %>% terra::plot()
#'  # there are tree crowns
#'  cloud2trees_ans$crowns_sf %>% dplyr::glimpse()
#'  cloud2trees_ans$crowns_sf %>% ggplot2::ggplot() + ggplot2::geom_sf(mapping = ggplot2::aes(fill = tree_height_m))
#'  # there are tree top points
#'  cloud2trees_ans$treetops_sf %>% dplyr::glimpse()
#'  cloud2trees_ans$treetops_sf %>% ggplot2::ggplot() + ggplot2::geom_sf(mapping = ggplot2::aes(color = tree_height_m))
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
  # message
    # start time
    xx1_cloud2raster <- Sys.time()
    message(
      "starting cloud2raster() step at ..."
      , xx1_cloud2raster
    )
  # do it
  cloud2raster_ans <- cloud2raster(
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
  # message
    # start time
    xx2_raster2trees <- Sys.time()
    message(
      "starting raster2trees() step at ..."
      , xx2_raster2trees
    )
  # do it
  raster2trees_ans <- raster2trees(
    chm_rast = cloud2raster_ans$chm_rast
    , outfolder = cloud2raster_ans$create_project_structure_ans$temp_dir
    , ws = ws
    , min_height = min_height
    , tempdir = cloud2raster_ans$create_project_structure_ans$temp_dir
  )

  if(nrow(raster2trees_ans)==0){
    stop(paste0(
      "No trees detected from point cloud data. Try adjusting some of these function parameters:"
      , "\nchm_res_m, min_height, max_height, ws"
    ))
  }

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
  # cloud2trees::trees_competition()
  ####################################################################
  # start time
  xx3_trees_competition <- Sys.time()
  if(estimate_tree_competition==T){
    # message
    message(
      "starting trees_competition() step at ..."
      , xx3_trees_competition
    )
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
  }else{
    # empty data
    trees_competition_ans <- dplyr::tibble(
      treeID = character(0)
      , comp_trees_per_ha = numeric(0)
      , comp_relative_tree_height = numeric(0)
      , comp_dist_to_nearest_m = numeric(0)
    )
  }

  ####################################################################
  # cloud2trees::treels_stem_dbh()
  ####################################################################
  # start time
  xx4_treels_stem_dbh <- Sys.time()
  if(estimate_dbh_from_cloud==T){
    # message
    message(
      "starting treels_stem_dbh() step at ..."
      , xx4_treels_stem_dbh
    )
    # do it
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
  # cloud2trees::trees_dbh()
  ####################################################################
  # start time
  xx5_trees_dbh <- Sys.time()
  if(estimate_dbh_from_cloud==T){
    # message
    message(
      "starting trees_dbh() step at ..."
      , xx5_trees_dbh
    )
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
    # message
    message(
      "starting trees_dbh() step at ..."
      , xx5_trees_dbh
    )
    # trees_dbh
    trees_dbh_ans <- trees_dbh(
      tree_list = raster2trees_ans
      , study_boundary = cloud2raster_ans$chunk_las_catalog_ans$las_ctg@data$geometry
      , dbh_model = dbh_model
      , treels_dbh_locations = NA
      , input_treemap_dir = cloud2raster_ans$create_project_structure_ans$input_treemap_dir
      , outfolder = cloud2raster_ans$create_project_structure_ans$delivery_dir
    )
  }else{
    # empty data
    trees_dbh_ans <- dplyr::tibble(
      treeID = character(0)
      , fia_est_dbh_cm = as.numeric(0)
      , fia_est_dbh_cm_lower = as.numeric(0)
      , fia_est_dbh_cm_upper = as.numeric(0)
      , dbh_cm = as.numeric(0)
      , is_training_data = as.logical(0)
      , dbh_m = as.numeric(0)
      , radius_m = as.numeric(0)
      , basal_area_m2 = as.numeric(0)
      , basal_area_ft2 = as.numeric(0)
      , ptcld_extracted_dbh_cm = as.numeric(0)
      , ptcld_predicted_dbh_cm = as.numeric(0)
    )
  }

  # trees_dbh_ans %>% class()
  # trees_dbh_ans %>% dplyr::glimpse()
  # trees_dbh_ans %>% dplyr::select(tidyselect::contains("dbh_cm")) %>% dplyr::glimpse()
  # trees_dbh_ans %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=dbh_cm))
  # trees_dbh_ans %>% ggplot2::ggplot(ggplot2::aes(x = tree_height_m, y = dbh_cm)) + ggplot2::geom_point()

  ####################################################################
  # write data
  ####################################################################
    # message
      # start time
      xx6_return <- Sys.time()
      message(
        "started writing final data at ..."
        , xx6_return
      )
    # do it
    # get names from dbh data
    dbh_names_temp <- c(
        "treeID"
        , get_list_diff(
          names(trees_dbh_ans %>% sf::st_drop_geometry())
          , names(raster2trees_ans %>% sf::st_drop_geometry())
        )
      )
    # get names from competition data
    comp_names_temp <- c(
        "treeID"
        , get_list_diff(
          names(trees_competition_ans %>% sf::st_drop_geometry())
          , names(raster2trees_ans %>% sf::st_drop_geometry())
        )
      )
    # join all data together for final return data
    crowns_sf_with_dbh <- raster2trees_ans %>%
      # join dbh data
      dplyr::left_join(
        trees_dbh_ans %>%
          sf::st_drop_geometry() %>%
          dplyr::select(dplyr::all_of(dbh_names_temp))
        , by = "treeID"
      ) %>%
      # join competition data
      dplyr::left_join(
        trees_competition_ans %>%
          sf::st_drop_geometry() %>%
          dplyr::select(dplyr::all_of(comp_names_temp))
        , by = "treeID"
      )

    ### write the data to the disk
    if(nrow(crowns_sf_with_dbh)>250e3){
      # split up the detected crowns
      crowns_sf_with_dbh <- crowns_sf_with_dbh %>%
        dplyr::arrange(as.numeric(tree_x),as.numeric(tree_y)) %>%
        # groups of 250k
        dplyr::mutate(grp = ceiling(dplyr::row_number()/250e3))

      write_fnl_temp <- crowns_sf_with_dbh$grp %>%
        unique() %>%
        purrr::map(function(x){
          ### write the data to the disk
          # crown vector polygons
          sf::st_write(
            crowns_sf_with_dbh %>%
              dplyr::filter(grp == x) %>%
              dplyr::select(-c(grp))
            , paste0(cloud2raster_ans$create_project_structure_ans$delivery_dir, "/final_detected_crowns_",x,".gpkg")
            , append = FALSE
            , quiet = TRUE
          )
          # tree top vector points
          sf::st_write(
            # get tree points
            crowns_sf_with_dbh %>%
              dplyr::filter(grp == x) %>%
              dplyr::select(-c(grp)) %>%
              sf::st_drop_geometry() %>%
              sf::st_as_sf(coords = c("tree_x", "tree_y"), crs = sf::st_crs(crowns_sf_with_dbh))
            , paste0(cloud2raster_ans$create_project_structure_ans$delivery_dir, "/final_detected_tree_tops_",x,".gpkg")
            , append = FALSE
            , quiet = TRUE
          )
          return(
            dplyr::tibble(
              crowns_file = paste0(cloud2raster_ans$create_project_structure_ans$delivery_dir, "/final_detected_crowns_",x,".gpkg")
              , trees_file = paste0(cloud2raster_ans$create_project_structure_ans$delivery_dir, "/final_detected_tree_tops_",x,".gpkg")
            )
          )
        }) %>%
        dplyr::bind_rows()
    }else{
        # crown vector polygons
        sf::st_write(
          crowns_sf_with_dbh
          , paste0(cloud2raster_ans$create_project_structure_ans$delivery_dir, "/final_detected_crowns.gpkg")
          , append = FALSE
          , quiet = TRUE
        )
        # tree top vector points
        sf::st_write(
          # get tree points
          crowns_sf_with_dbh %>%
            sf::st_drop_geometry() %>%
            sf::st_as_sf(coords = c("tree_x", "tree_y"), crs = sf::st_crs(crowns_sf_with_dbh))
          , paste0(cloud2raster_ans$create_project_structure_ans$delivery_dir, "/final_detected_tree_tops.gpkg")
          , append = FALSE
          , quiet = TRUE
        )
    }

    # tree top points for return
    treetops_sf_with_dbh <- crowns_sf_with_dbh %>%
      sf::st_drop_geometry() %>%
      sf::st_as_sf(coords = c("tree_x", "tree_y"), crs = sf::st_crs(crowns_sf_with_dbh))

    # remove temp files
    if(keep_intrmdt==F){
      unlink(cloud2raster_ans$create_project_structure_ans$temp_dir, recursive = T)
    }

    # write settings and timer data
    xx7_fin <- Sys.time()
    # data
      return_df <-
          # data from las_ctg
          cloud2raster_ans$chunk_las_catalog_ans$las_ctg@data %>%
          sf::st_set_geometry("geometry") %>%
          dplyr::summarise(
            geometry = sf::st_union(geometry)
            , number_of_points = sum(Number.of.point.records, na.rm = T)
          ) %>%
          dplyr::mutate(
            las_area_m2 = sf::st_area(geometry) %>% as.numeric()
          ) %>%
          sf::st_drop_geometry() %>%
          dplyr::mutate(
            timer_cloud2raster_mins = difftime(xx2_raster2trees, xx1_cloud2raster, units = c("mins")) %>%
              as.numeric()
            , timer_raster2trees_mins = difftime(xx3_trees_competition, xx2_raster2trees, units = c("mins")) %>%
              as.numeric()
            , timer_trees_competition_mins = difftime(xx4_treels_stem_dbh, xx3_trees_competition, units = c("mins")) %>%
              as.numeric()
            , timer_treels_stem_dbh_mins = difftime(xx5_trees_dbh, xx4_treels_stem_dbh, units = c("mins")) %>%
              as.numeric()
            , timer_trees_dbh_mins = difftime(xx6_return, xx5_trees_dbh, units = c("mins")) %>%
              as.numeric()
            , timer_write_data_mins = difftime(xx7_fin, xx6_return, units = c("mins")) %>%
              as.numeric()
            , timer_total_time_mins = difftime(xx7_fin, xx1_cloud2raster, units = c("mins")) %>%
              as.numeric()
            # settings
            , sttng_input_las_dir = input_las_dir
            , sttng_accuracy_level = accuracy_level
            , sttng_max_ctg_pts = max_ctg_pts
            , sttng_max_area_m2 = max_area_m2
            , sttng_dtm_res_m = dtm_res_m
            , sttng_chm_res_m = chm_res_m
            , sttng_min_height = min_height
            , sttng_max_height = max_height
            , sttng_ws = as.character(deparse(ws)) %>% paste(collapse = "") %>% stringr::str_trim()
            , sttng_estimate_tree_dbh = estimate_tree_dbh
            , sttng_max_dbh = max_dbh
            , sttng_dbh_model = dbh_model
            , sttng_estimate_dbh_from_cloud = estimate_dbh_from_cloud
            , sttng_estimate_tree_competition = estimate_tree_competition
            , sttng_competition_buffer_m = competition_buffer_m
            , sttng_search_dist_max = search_dist_max
          )

      # write
      write.csv(
        return_df
        , paste0(cloud2raster_ans$create_project_structure_ans$delivery_dir, "/processed_tracking_data.csv")
        , row.names = F
      )
  ####################################################################
  # return
  ####################################################################
    # message
    message(
      "cloud2trees() total time was "
      , round(as.numeric(difftime(xx7_fin, xx1_cloud2raster, units = c("mins"))),2)
      , " minutes to process "
      , scales::comma(sum(cloud2raster_ans$chunk_las_catalog_ans$las_ctg@data$Number.of.point.records))
      , " points over an area of "
      , scales::comma(as.numeric(cloud2raster_ans$chunk_las_catalog_ans$las_ctg@data$geometry %>% sf::st_union() %>% sf::st_area())/10000,accuracy = 0.01)
      , " hectares"
    )
    # return
    return(list(
      crowns_sf = crowns_sf_with_dbh
      , treetops_sf = treetops_sf_with_dbh
      , dtm_rast = cloud2raster_ans$dtm_rast
      , chm_rast = cloud2raster_ans$chm_rast
    ))
}

####################################################################
## intermediate fn to get list difference
####################################################################
get_list_diff <- function(x, y) {
  if(inherits(x, "character") & inherits(y, "character")){
      d <- x[!(x %in% y)]
      d <- unique(d)
  }else{
    stop("must provide character list in x and y")
  }
  return(d)
}
