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
#' * Extract tree forest type group using [trees_type()] (*if set to TRUE*)
#' * Extract tree CBH values from the normalized point cloud and estimate missing values using [trees_cbh()] (*if set to TRUE*)
#' * Estimate tree biomass (or crown biomass) using [trees_biomass()] (*if method is denoted*)
#'
#' See the documentation for each individual function called for more details.
#'
#' @inheritParams cloud2raster
#' @inheritParams find_ext_data
#' @param ws numeric or function. Length or diameter of the moving window used to detect the local
#' maxima in the units of the input data (usually meters). If it is numeric a fixed window size is used.
#' If it is a function, the function determines the size of the window at any given location on the canopy.
#' By default function takes the height of a given pixel as its only argument and return the
#' desired size of the search window when centered on that pixel.
#' @param estimate_tree_dbh logical. Should tree DBH be estimated? See [trees_dbh()].
#' @param max_dbh numeric. Set the largest tree diameter (m) expected in the point cloud
#' @param dbh_model `r lifecycle::badge("deprecated")` Use the `dbh_model_regional` or `dbh_model_local` argument instead.
#' @param dbh_model_regional string. Set the model to use for regional dbh-height allometry based on FIA tree measurements.
#' Can be "cr" for the Chapman-Richards formula (default) or "power" for power function
#' @param dbh_model_local string. Set the model to use for local dbh-height allometry based on provided DBH training data in `treels_dbh_locations`.
#' Can be "rf" for random forest or "lin" for linear
#' @param estimate_dbh_from_cloud logical. Should DBH be estimated from the point cloud? See [treels_stem_dbh()]. Setting to `TRUE` may significantly increase processing time.
#' @param estimate_tree_competition logical. Should tree competition metrics be calculated? See [trees_competition()]. Setting to `TRUE` may slightly increase processing time.
#' @param competition_buffer_m number. Set buffer around tree (m) to calculate competition metrics
#' @param search_dist_max `r lifecycle::badge("deprecated")` Use the `competition_max_search_dist_m` argument instead.
#' @param competition_max_search_dist_m number. Maximum search distance (m) to nearest tree for competition. Larger search distances will increase processing time and possibly result in memory issues.
#' If no competition trees are found within this distance, the return column `comp_dist_to_nearest_m` = `competition_max_search_dist_m` parameter.
#' @param estimate_tree_type logical. Should tree forest type be estimated? See [trees_type()].
#' @param type_max_search_dist_m number. Maximum search distance (m) to obtain forest type group data for trees that overlap with non-forest data in the original Wilson (2023) data.
#' Larger search distances will increase processing time and possibly result in memory issues.
#' @param estimate_tree_hmd logical. Should tree height of the maximum crown diameter (HMD) be estimated? See [trees_hmd()].
#' @param hmd_tree_sample_n,hmd_tree_sample_prop numeric. Provide either `tree_sample_n`, the number of trees, or `tree_sample_prop`, the
#'   proportion of the trees to attempt to extract a HMD from the point cloud for.
#'   If neither are supplied, `tree_sample_n = 777` will be used. If both are supplied, `tree_sample_n` will be used.
#'   Increasing `tree_sample_prop` toward one (1) will increase the processing time, perhaps significantly depending on the number of trees in the `trees_poly` data.
#'   The maximum number of trees to extract tree HMD using `cloud2trees()` is 20,000.
#'   Try `trees_hmd()` with outputs from `cloud2trees()` if you want to attempt to extract HMD for >20,000 trees.
#' @param hmd_estimate_missing_hmd logical. It is not likely that HMD will be extracted successfully from every tree.
#'   Should the missing HMD values be estimated using the tree height and location information based on trees for which HMD is successfully extracted?
#' @param estimate_biomass_method character. To estimate tree biomass or tree (or crown biomass) enter one or a list of multiple biomass methods. See [trees_biomass()].
#'   Leave as blank (i.e. `NA`) to skip biomass estimation.
#' @param biomass_max_crown_kg_per_m3 numeric. the maximum CBD of the tree crown in kilograms per cubic meter.
#' Values above this limit will be set at the median value for the area using only stands that have CBD values lower than this limit.
#' The default value of 2 kilograms per cubic meter was based on [Mell et al. (2009)](https://doi.org/10.1016/j.combustflame.2009.06.015)
#' who found the dry bulk density of the tree crown was 2.6 kilograms per cubed meter
#' using Douglas-fir trees grown on Christmas tree farms.
#' Set this parameter to a large value (e.g. 1e10) or NULL to avoid limiting tree crown CBD.
#' @param estimate_tree_cbh logical. Should tree DBH be estimated? See [trees_cbh()].
#'   Make sure to set `cbh_estimate_missing_cbh = TRUE` if you want to obtain CBH values for cases when CBH cannot be extracted from the point cloud.
#' @param cbh_tree_sample_n,cbh_tree_sample_prop numeric. Provide either `tree_sample_n`, the number of trees, or `tree_sample_prop`, the
#'   proportion of the trees to attempt to extract a CBH from the point cloud for.
#'   If neither are supplied, `tree_sample_n = 333` will be used. If both are supplied, `tree_sample_n` will be used.
#'   Increasing `tree_sample_prop` toward one (1) will increase the processing time, perhaps significantly depending on the number of trees in the `trees_poly` data.
#'   The maximum number of trees to extract tree CBH using `cloud2trees()` is 20,000.
#'   Try `trees_cbh()` with outputs from `cloud2trees()` if you want to attempt to extract CBH for >20,000 trees.
#' @param cbh_which_cbh character. One of: "lowest"; "highest"; or "max_lad". See Viedma et al. (2024) reference.
#' * "lowest" - Height of the CBH of the segmented tree based on the last distance found in its profile
#' * "highest" - Height of the CBH of the segmented tree based on the maximum distance found in its profile
#' * "max_lad" - Height of the CBH of the segmented tree based on the maximum LAD percentage
#' @param cbh_estimate_missing_cbh logical. even if the `cbh_tree_sample_prop` parameter is set to "1", it is not likely that CBH will be extracted successfully from every tree.
#'   Should the missing CBH values be estimated using the tree height and location information based on trees for which CBH is successfully extracted?
#' @param cbh_min_vhp_n numeric. the minimum number of vertical height profiles (VHPs) needed to estimate a CBH.
#' @param cbh_voxel_grain_size_m numeric. horizontal resolution (suggested 1 meter for lad profiles and 10 meters for LAI maps). See `grain.size` in [leafR::lad.voxels()]
#' @param cbh_dist_btwn_bins_m numeric. value for the actual height bin step (in meters). See `step` in [LadderFuelsR::get_gaps_fbhs()]
#' @param cbh_min_fuel_layer_ht_m numeric. value for the actual minimum base height (in meters). See `min_height` in [LadderFuelsR::get_gaps_fbhs()]
#' @param cbh_lad_pct_gap numeric. value of the percentile threshold used to identify gaps (default percentile 25th). See `perc_gap` in [LadderFuelsR::get_gaps_fbhs()]
#' @param cbh_lad_pct_base numeric. value of the percentile threshold used to identify fuels layers base height (default percentile 25th). See `perc_base` in [LadderFuelsR::get_gaps_fbhs()]
#' @param cbh_num_jump_steps numeric. value for the number of height bin steps that can be jumped to reshape fuels layers. See `number_steps` in [LadderFuelsR::get_real_fbh()]
#' @param cbh_min_lad_pct numeric. value for the minimum required LAD percentage in a fuel layer. See `threshold` in [LadderFuelsR::get_layers_lad()]
#' @param cbh_frst_layer_min_ht_m numeric. value for the depth height of the first fuel layer. If the first fuel layer has the maximum LAD and its depth is greater than the indicated value, then this fuel layer is considered as the CBH of the tree. On the contrary, if its depth is <= the value, the CBH with maximum LAD will be the second fuel layer, although it has not the maximum LAD. See `hdepth1_height` in [LadderFuelsR::get_cbh_metrics()]
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
#'  cloud2trees_ans$crowns_sf %>% ggplot2::ggplot() +
#'   ggplot2::geom_sf(mapping = ggplot2::aes(fill = tree_height_m))
#'  # there are tree top points
#'  cloud2trees_ans$treetops_sf %>% dplyr::glimpse()
#'  cloud2trees_ans$treetops_sf %>% ggplot2::ggplot() +
#'   ggplot2::geom_sf(mapping = ggplot2::aes(color = tree_height_m))
#'  }
#' @export
#'
cloud2trees <- function(
  output_dir
  , input_las_dir
  , input_treemap_dir = NULL
  , input_foresttype_dir = NULL
  , input_landfire_dir = NULL
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
  , ws = itd_ws_functions()[["log_fn"]]
  , estimate_tree_dbh = FALSE
  , max_dbh = 2
  , dbh_model_regional = "cr"
  , dbh_model_local = "lin"
  , estimate_dbh_from_cloud = FALSE
  , estimate_tree_competition = FALSE
  , competition_buffer_m = 5
  , search_dist_max
  , competition_max_search_dist_m = 10
  , estimate_tree_type = FALSE
  , type_max_search_dist_m = 1000
  , estimate_tree_hmd = FALSE
  , hmd_tree_sample_n = NA
  , hmd_tree_sample_prop = NA
  , hmd_estimate_missing_hmd = FALSE
  , estimate_biomass_method = NA
  , biomass_max_crown_kg_per_m3 = 2
  , estimate_tree_cbh = FALSE
  , cbh_tree_sample_n = NA
  , cbh_tree_sample_prop = NA
  , cbh_which_cbh = "lowest"
  , cbh_estimate_missing_cbh = FALSE
  , cbh_min_vhp_n = 3
  , cbh_voxel_grain_size_m = 1
  , cbh_dist_btwn_bins_m = 1
  , cbh_min_fuel_layer_ht_m = 1
  , cbh_lad_pct_gap = 25
  , cbh_lad_pct_base = 25
  , cbh_num_jump_steps = 1
  , cbh_min_lad_pct = 10
  , cbh_frst_layer_min_ht_m = 1
  , overwrite = TRUE
){
  # check ws
  ws <- check_numeric_returning_function(ws)

  ####################################################################
  # check deprecated parameters
  ####################################################################
    calls <- names(sapply(match.call(), deparse))[-1]
    if(any("search_dist_max" %in% calls)) {
        stop(
          "`search_dist_max` deprecated. Use the `competition_max_search_dist_m` argument instead."
        )
    }
    if(any("dbh_model" %in% calls)) {
        stop(
          "`dbh_model` deprecated. Use the `dbh_model_regional` or `dbh_model_local` argument instead."
        )
    }
  ####################################################################
  # check biomass method so we can throw error before we kick off
  ####################################################################
    if(
      any(!is.na(estimate_biomass_method)) && any(!is.null(estimate_biomass_method))
      && is.character(estimate_biomass_method)
    ){
      which_biomass_methods <- check_biomass_method(estimate_biomass_method)
    }else{
      which_biomass_methods <- NULL
    }

    #### we now need to ensure we get DBH and CBH if we want biomass
    if(!is.null(which_biomass_methods)){
      estimate_tree_cbh <- T
      cbh_estimate_missing_cbh <- T
      estimate_tree_dbh <- T
    }
    if(any(stringr::str_equal(which_biomass_methods, "cruz"))){
      estimate_tree_type <- T
    }
  ####################################################################
  # check external data
  ####################################################################
    # find external data
    find_ext_data_ans <- find_ext_data(
      input_treemap_dir = input_treemap_dir
      , input_foresttype_dir = input_foresttype_dir
      , input_landfire_dir = input_landfire_dir
    )
    # if can't find external treemap data
    if(
      (estimate_dbh_from_cloud == T || estimate_tree_dbh == T) &&
      is.null(find_ext_data_ans$treemap_dir)
    ){
      stop(paste0(
        "Treemap data has not been downloaded to package contents. Use `get_treemap()` first."
        , "\nIf you supplied a value to the `input_treemap_dir` parameter check that directory for data."
      ))
    }
    # if can't find external foresttype data
    if(estimate_tree_type == T && is.null(find_ext_data_ans$foresttype_dir)){
      stop(paste0(
        "Forest Type Group data has not been downloaded to package contents. Use `get_foresttype()` first."
        , "\nIf you supplied a value to the `input_foresttype_dir` parameter check that directory for data."
      ))
    }
    # if can't find external landfire data
    if(
      any(stringr::str_equal(which_biomass_methods, "landfire")) &&
      is.null(find_ext_data_ans$landfire_dir)
    ){
      stop(paste0(
        "LANDFIRE data has not been downloaded to package contents. Use `get_landfire()` first."
        , "\nIf you supplied a value to the `input_landfire_dir` parameter check that directory for data."
      ))
    }

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
  # update intermediate keep
  updt_keep_intrmdt <- dplyr::case_when(
    keep_intrmdt==F & (estimate_tree_cbh==T | estimate_tree_hmd==T | estimate_dbh_from_cloud==T) ~ T
    , T ~ keep_intrmdt
  )
  # do it
  cloud2raster_ans <- cloud2raster(
      output_dir = output_dir
      , input_las_dir = input_las_dir
      , input_treemap_dir = find_ext_data_ans$treemap_dir
      , input_foresttype_dir = find_ext_data_ans$foresttype_dir
      , accuracy_level = accuracy_level
      , max_ctg_pts = max_ctg_pts
      , max_area_m2 = max_area_m2
      , transform = transform
      , new_crs = new_crs
      , old_crs = old_crs
      , keep_intrmdt = updt_keep_intrmdt
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
    , outfolder = cloud2raster_ans$create_project_structure_ans$delivery_dir # write crown polygons and tree top points
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
  err_trees_competition <- NULL
  # empty data
    trees_competition_ans_temp <- dplyr::tibble(
      treeID = character(0)
      , comp_trees_per_ha = numeric(0)
      , comp_relative_tree_height = numeric(0)
      , comp_dist_to_nearest_m = numeric(0)
    )
  if(estimate_tree_competition==T){
    # message
    message(
      "starting trees_competition() step at ..."
      , xx3_trees_competition
    )
    # trees_competition
    # safe it
    safe_trees_competition <- purrr::safely(trees_competition)
    trees_competition_ans <- safe_trees_competition(
      tree_list = raster2trees_ans
      , competition_buffer_m = competition_buffer_m
      , study_boundary = cloud2raster_ans$chunk_las_catalog_ans$las_ctg@data$geometry
      , search_dist_max = as.numeric(competition_max_search_dist_m)
    )
    # handle error
    if(is.null(trees_competition_ans$error)){ # no error
      # just get the result
      trees_competition_ans <- trees_competition_ans$result
    }else{
      # error
      err_trees_competition <- trees_competition_ans$error
      # empty data
      trees_competition_ans <- trees_competition_ans_temp
    }

    # trees_competition_ans %>% class()
    # trees_competition_ans %>% dplyr::glimpse()
    # trees_competition_ans %>% dplyr::select(tidyselect::starts_with("comp_")) %>% dplyr::glimpse()
    # trees_competition_ans %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=tree_height_m))
    # trees_competition_ans %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=comp_dist_to_nearest_m))
  }else{
    # empty data
    trees_competition_ans <- trees_competition_ans_temp
  }

  ####################################################################
  # cloud2trees::treels_stem_dbh()
  ####################################################################
  # start time
  xx4_treels_stem_dbh <- Sys.time()
  err_treels_stem_dbh <- NULL
  if(estimate_dbh_from_cloud==T){
    # message
    message(
      "starting treels_stem_dbh() step at ..."
      , xx4_treels_stem_dbh
    )
    # do it
    # treels stuff
    safe_treels_stem_dbh <- purrr::safely(treels_stem_dbh)
    treels_dbh_locations <- safe_treels_stem_dbh(
      folder = cloud2raster_ans$normalize_flist
      , outfolder = cloud2raster_ans$create_project_structure_ans$treels_dir
      , min_height = min_height
      , max_dbh = max_dbh
      , chunk_these = !cloud2raster_ans$chunk_las_catalog_ans$is_chunked_grid
    )

    # handle error
    if(is.null(treels_dbh_locations$error)){ # no error
      # just get the result
      treels_dbh_locations <- treels_dbh_locations$result
    }else{
      # error
      err_treels_stem_dbh <- treels_dbh_locations$error
      # empty data
      treels_dbh_locations <- NA
    }
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
  err_trees_dbh <- NULL
  # empty data
    trees_dbh_ans_temp <- dplyr::tibble(
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
  if(estimate_dbh_from_cloud==T){
    # message
    message(
      "starting trees_dbh() step at ..."
      , xx5_trees_dbh
    )
    # trees_dbh
    safe_trees_dbh <- purrr::safely(trees_dbh)
    trees_dbh_ans <- safe_trees_dbh(
      tree_list = raster2trees_ans
      , study_boundary = cloud2raster_ans$chunk_las_catalog_ans$las_ctg@data$geometry
      , dbh_model_regional = dbh_model_regional
      , dbh_model_local = dbh_model_local
      , treels_dbh_locations = treels_dbh_locations
      , input_treemap_dir = cloud2raster_ans$create_project_structure_ans$input_treemap_dir
      , outfolder = cloud2raster_ans$create_project_structure_ans$delivery_dir
    )
    # handle error
    if(is.null(trees_dbh_ans$error)){ # no error
      # just get the result
      trees_dbh_ans <- trees_dbh_ans$result
    }else{
      # error
      err_trees_dbh <- trees_dbh_ans$error
      # empty data
      trees_dbh_ans <- trees_dbh_ans_temp
    }
  }else if(estimate_tree_dbh==T){
    # message
    message(
      "starting trees_dbh() step at ..."
      , xx5_trees_dbh
    )
    # trees_dbh
    safe_trees_dbh <- purrr::safely(trees_dbh)
    trees_dbh_ans <- safe_trees_dbh(
      tree_list = raster2trees_ans
      , study_boundary = cloud2raster_ans$chunk_las_catalog_ans$las_ctg@data$geometry
      , dbh_model_regional = dbh_model_regional
      , dbh_model_local = dbh_model_local
      , treels_dbh_locations = NA
      , input_treemap_dir = cloud2raster_ans$create_project_structure_ans$input_treemap_dir
      , outfolder = cloud2raster_ans$create_project_structure_ans$delivery_dir
    )
    # handle error
    if(is.null(trees_dbh_ans$error)){ # no error
      # just get the result
      trees_dbh_ans <- trees_dbh_ans$result
    }else{
      # error
      err_trees_dbh <- trees_dbh_ans$error
      # empty data
      trees_dbh_ans <- trees_dbh_ans_temp
    }
  }else{
    # empty data
    trees_dbh_ans <- trees_dbh_ans_temp
  }

  # trees_dbh_ans %>% class()
  # trees_dbh_ans %>% dplyr::glimpse()
  # trees_dbh_ans %>% dplyr::select(tidyselect::contains("dbh_cm")) %>% dplyr::glimpse()
  # trees_dbh_ans %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=dbh_cm))
  # trees_dbh_ans %>% ggplot2::ggplot(ggplot2::aes(x = tree_height_m, y = dbh_cm)) + ggplot2::geom_point()

  ####################################################################
  # cloud2trees::trees_cbh()
  ####################################################################
  # start time
  xx6_trees_cbh <- Sys.time()
  err_trees_cbh <- NULL
  err_trees_cbh_sample <- F
  # empty data
    trees_cbh_ans_temp <- dplyr::tibble(
      treeID = character(0)
      , tree_cbh_m = as.numeric(0)
      , is_training_cbh = as.logical(0)
    )
  if(estimate_tree_cbh==T){
    # message
    message(
      "starting trees_cbh() step at ..."
      , xx6_trees_cbh
    )
    ### limit cbh sample to 20,000
    ### if one wants more, they can run trees_cbh in standalone
    if(
      dplyr::coalesce(cbh_tree_sample_n,0)>20000
      || dplyr::coalesce(cbh_tree_sample_prop,0)*nrow(raster2trees_ans)>20000
    ){
      cbh_tree_sample_n <- 20000
      err_trees_cbh_sample <- T
      message(paste0(
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        , "\n"
        , "The maximum number of trees to extract tree CBH using cloud2trees() is 20,000"
        , "\n..............try trees_cbh() if want to attempt to extract CBH for >20,000 trees"
        , "\nusing `final_detected_crowns.gpkg` and `norm_las` in the `point_cloud_processing_delivery` directory"
      ))
    }
    # trees_cbh
    safe_trees_cbh <- purrr::safely(trees_cbh)
    trees_cbh_ans <- safe_trees_cbh(
      trees_poly = raster2trees_ans
      , norm_las = cloud2raster_ans$create_project_structure_ans$las_normalize_dir
      , tree_sample_n = cbh_tree_sample_n
      , tree_sample_prop = cbh_tree_sample_prop
      , which_cbh = cbh_which_cbh
      , estimate_missing_cbh = cbh_estimate_missing_cbh
      , min_vhp_n = cbh_min_vhp_n
      , voxel_grain_size_m = cbh_voxel_grain_size_m
      , dist_btwn_bins_m = cbh_dist_btwn_bins_m
      , min_fuel_layer_ht_m = cbh_min_fuel_layer_ht_m
      , lad_pct_gap = cbh_lad_pct_gap
      , lad_pct_base = cbh_lad_pct_base
      , num_jump_steps = cbh_num_jump_steps
      , min_lad_pct = cbh_min_lad_pct
      , frst_layer_min_ht_m = cbh_frst_layer_min_ht_m
      , force_same_crs = T
      , outfolder = cloud2raster_ans$create_project_structure_ans$delivery_dir
    )
    # handle error
    if(is.null(trees_cbh_ans$error)){ # no error
      # just get the result
      trees_cbh_ans <- trees_cbh_ans$result
    }else{
      # error
      err_trees_cbh <- trees_cbh_ans$error
      # empty data
      trees_cbh_ans <- trees_cbh_ans_temp
    }
  }else{
    # empty data
    trees_cbh_ans <- trees_cbh_ans_temp
  }

  ####################################################################
  # cloud2trees::trees_type()
  ####################################################################
  # start time
  xx7_trees_type <- Sys.time()
  err_trees_type <- NULL
  # empty data
    trees_type_ans_temp <- dplyr::tibble(
      treeID = character(0)
      , forest_type_group_code = character(0)
      , forest_type_group = character(0)
      , hardwood_softwood = character(0)
    )
    trees_type_rast <- NULL
  if(estimate_tree_type==T){
    # message
    message(
      "starting trees_type() step at ..."
      , xx7_trees_type
    )
    # trees_type
    safe_trees_type <- purrr::safely(trees_type)
    trees_type_ans <- safe_trees_type(
      tree_list = raster2trees_ans
      , study_boundary = cloud2raster_ans$chunk_las_catalog_ans$las_ctg@data$geometry
      , input_foresttype_dir = cloud2raster_ans$create_project_structure_ans$input_foresttype_dir
      , max_search_dist_m = as.numeric(type_max_search_dist_m)
    )
    # handle error
    if(is.null(trees_type_ans$error)){ # no error
      # just get the result
      trees_type_rast <- trees_type_ans$result$foresttype_rast # $foresttype_rast unique to this function
      trees_type_ans <- trees_type_ans$result$tree_list  # $tree_list unique to this function
    }else{
      # error
      err_trees_type <- trees_type_ans$error
      # empty data
      trees_type_ans <- trees_type_ans_temp
    }
  }else{
    # empty data
    trees_type_ans <- trees_type_ans_temp
  }

  ####################################################################
  # cloud2trees::trees_hmd()
  ####################################################################
  # start time
  xx8_trees_hmd <- Sys.time()
  err_trees_hmd <- NULL
  err_trees_hmd_sample <- F
  # empty data
    trees_hmd_ans_temp <- dplyr::tibble(
      treeID = character(0)
      , max_crown_diam_height_m = as.numeric(0)
      , is_training_hmd = as.logical(0)
    )
  if(estimate_tree_hmd==T){
    # message
    message(
      "starting trees_hmd() step at ..."
      , xx8_trees_hmd
    )
    ### limit hmd sample to 20,000
    ### if one wants more, they can run trees_hmd in standalone
    if(
      dplyr::coalesce(hmd_tree_sample_n,0)>20000
      || dplyr::coalesce(hmd_tree_sample_prop,0)*nrow(raster2trees_ans)>20000
    ){
      hmd_tree_sample_n <- 20000
      err_trees_hmd_sample <- T
      message(paste0(
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        , "\n"
        , "The maximum number of trees to extract tree HMD using cloud2trees() is 20,000"
        , "\n..............try trees_hmd() if want to attempt to extract HMD for >20,000 trees"
        , "\nusing `final_detected_crowns.gpkg` and `norm_las` in the `point_cloud_processing_delivery` directory"
      ))
    }
    # trees_hmd
    safe_trees_hmd <- purrr::safely(trees_hmd)
    trees_hmd_ans <- safe_trees_hmd(
      trees_poly = raster2trees_ans
      , norm_las = cloud2raster_ans$create_project_structure_ans$las_normalize_dir
      , tree_sample_n = hmd_tree_sample_n
      , tree_sample_prop = hmd_tree_sample_prop
      , estimate_missing_hmd = hmd_estimate_missing_hmd
      , force_same_crs = T
      , outfolder = cloud2raster_ans$create_project_structure_ans$delivery_dir
    )
    # handle error
    if(is.null(trees_hmd_ans$error)){ # no error
      # just get the result
      trees_hmd_ans <- trees_hmd_ans$result
    }else{
      # error
      err_trees_hmd <- trees_hmd_ans$error
      # empty data
      trees_hmd_ans <- trees_hmd_ans_temp
    }
  }else{
    # empty data
    trees_hmd_ans <- trees_hmd_ans_temp
  }

  ####################################################################
  # cloud2trees::trees_biomass()
  ####################################################################
  # start time
  xx9_trees_biomass <- Sys.time()
  err_trees_biomass <- NULL
  # empty data
    trees_biomass_ans_temp <- dplyr::tibble(
      treeID = character(0)
    )
    trees_biomass_cell_landfire <- NULL
    trees_biomass_cell_cruz <- NULL
  ### we need to make sure:
    ### 1. we want biomass; 2. we have dbh; 3. we have cbh
  if(
    !is.null(which_biomass_methods) &&
    is.null(err_trees_cbh) &&
    is.null(err_trees_dbh)
  ){
    # message
    message(
      "starting trees_biomass() step at ..."
      , xx9_trees_biomass
    )
    # check if already done trees_type() so we don't do it twice and
      # potentially get different answer
    if(
      is.null(err_trees_type) &&
      inherits(trees_type_ans, "data.frame")
    ){
      # check for forest_type exists and not all missing
      safe_check_df_cols_all_missing <- purrr::safely(check_df_cols_all_missing)
      chk_foresttype <- safe_check_df_cols_all_missing(
          trees_type_ans
          , col_names = c("forest_type_group_code", "forest_type_group")
          , all_numeric = F
        )
      # build join data
      if(is.null(chk_foresttype$error)){
        # columns added by trees_type()
        type_names_temp <- c(
          "treeID"
          , get_list_diff(
            names(trees_type_ans %>% sf::st_drop_geometry())
            , names(raster2trees_ans %>% sf::st_drop_geometry())
          )
        )
        # join data
        df_foresttype <- trees_type_ans %>%
          sf::st_drop_geometry() %>%
          dplyr::select(dplyr::all_of(type_names_temp))
      }else{
        # blank data
        df_foresttype <- dplyr::tibble(treeID = character(0))
      }

    }else{
      # blank data
      df_foresttype <- dplyr::tibble(treeID = character(0))
    }

    # trees_biomass
    safe_trees_biomass <- purrr::safely(trees_biomass)
    trees_biomass_ans <- safe_trees_biomass(
      #### !!!!!!!!!!!!!!!!!!!!!! NOTICE WE HAVE TO JOIN THE TREE LIST HERE
      tree_list = raster2trees_ans %>%
        # join foresttype
        dplyr::left_join(
          df_foresttype
          , by = "treeID"
        ) %>%
        # join dbh data
        dplyr::left_join(
          trees_dbh_ans %>%
            sf::st_drop_geometry() %>%
            dplyr::select(dplyr::all_of(c("treeID", "dbh_cm")))
          , by = "treeID"
        ) %>%
        # join cbh data
        dplyr::left_join(
          trees_cbh_ans %>%
            sf::st_drop_geometry() %>%
            dplyr::select(dplyr::all_of(c("treeID", "tree_cbh_m")))
          , by = "treeID"
        )
      , study_boundary = cloud2raster_ans$chunk_las_catalog_ans$las_ctg@data$geometry
      , method = which_biomass_methods
      , max_crown_kg_per_m3 = biomass_max_crown_kg_per_m3
    )
    # handle error
    if(is.null(trees_biomass_ans$error)){ # no error
      # stand_cell_data_* result
      trees_biomass_cell_landfire <- trees_biomass_ans$result$stand_cell_data_landfire # $stand_cell_data_landfire unique to this function
      trees_biomass_cell_cruz <- trees_biomass_ans$result$stand_cell_data_cruz # $stand_cell_data_landfire unique to this function
      # just get the result
      trees_biomass_ans <- trees_biomass_ans$result$tree_list %>% # $tree_list unique to this function
        dplyr::select( -dplyr::any_of(c(
          "hey_xxxxxxxxxx"
          , "dbh_m"
          , "dbh_cm"
          , "basal_area_m2"
          , "tree_cbh_m"
          , "forest_type_group_code"
          , "forest_type_group"
          , "hardwood_softwood"
        )))

    }else{
      # error
      err_trees_biomass <- trees_biomass_ans$error
      # empty data
      trees_biomass_ans <- trees_biomass_ans_temp
    }
  }else{
    # empty data
    trees_biomass_ans <- trees_biomass_ans_temp
    # there were dbh, cbh errors
    if(!is.null(which_biomass_methods)){
      err_trees_biomass <- "Error: could not execute trees_biomass() step...see DBH and/or CBH error"
    }
  }

  ####################################################################
  # write data
  ####################################################################
    # message
      # start time
      xx10_return <- Sys.time()
      message(
        "started writing final data at ..."
        , xx10_return
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
    # get names from cbh data
    cbh_names_temp <- c(
        "treeID"
        , get_list_diff(
          names(trees_cbh_ans %>% sf::st_drop_geometry())
          , names(raster2trees_ans %>% sf::st_drop_geometry())
        )
      )
    # get names from hmd data
    hmd_names_temp <- c(
        "treeID"
        , get_list_diff(
          names(trees_hmd_ans %>% sf::st_drop_geometry())
          , names(raster2trees_ans %>% sf::st_drop_geometry())
        )
      )
    # get names from type data
    type_names_temp <- c(
        "treeID"
        , get_list_diff(
          names(trees_type_ans %>% sf::st_drop_geometry())
          , names(raster2trees_ans %>% sf::st_drop_geometry())
        )
      )
    # get names from biomass data
    biomass_names_temp <- c(
        "treeID"
        , get_list_diff(
          names(trees_biomass_ans %>% sf::st_drop_geometry())
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
      # join cbh data
      dplyr::left_join(
        trees_cbh_ans %>%
          sf::st_drop_geometry() %>%
          dplyr::select(dplyr::all_of(cbh_names_temp))
        , by = "treeID"
      ) %>%
      # join type data
      dplyr::left_join(
        trees_type_ans %>%
          sf::st_drop_geometry() %>%
          dplyr::select(dplyr::all_of(type_names_temp))
        , by = "treeID"
      ) %>%
      # join competition data
      dplyr::left_join(
        trees_competition_ans %>%
          sf::st_drop_geometry() %>%
          dplyr::select(dplyr::all_of(comp_names_temp))
        , by = "treeID"
      ) %>%
      # join hmd data
      dplyr::left_join(
        trees_hmd_ans %>%
          sf::st_drop_geometry() %>%
          dplyr::select(dplyr::all_of(hmd_names_temp))
        , by = "treeID"
      ) %>%
      # join biomass data
      dplyr::left_join(
        trees_biomass_ans %>%
          sf::st_drop_geometry() %>%
          dplyr::select(dplyr::all_of(biomass_names_temp))
        , by = "treeID"
      )

    ### write crown polygons and tree top points
    ### input data should be from raster2trees and can include extra columns
    ### !!!! this will overwrite files already written above
    write_fnl <- write_raster2trees_ans(
      raster2trees_ans = crowns_sf_with_dbh
      , dir = cloud2raster_ans$create_project_structure_ans$delivery_dir
    )

    # tree top points for return
    treetops_sf_with_dbh <- crowns_sf_with_dbh %>%
      sf::st_drop_geometry() %>%
      sf::st_as_sf(coords = c("tree_x", "tree_y"), crs = sf::st_crs(crowns_sf_with_dbh))

    # remove temp files
    if(
      keep_intrmdt==F &&
      estimate_tree_cbh==F &&
      estimate_tree_hmd==F &&
      estimate_dbh_from_cloud==F &&
      length(list.files(cloud2raster_ans$create_project_structure_ans$reproj_dir))==0
    ){
      unlink(cloud2raster_ans$create_project_structure_ans$temp_dir, recursive = T)
    }else if(
      (estimate_tree_cbh==T || estimate_tree_hmd==T || estimate_dbh_from_cloud==T)
    ){
      # move the normalized laz files to delivery
      ### Create the directories
      to_dir <- file.path(cloud2raster_ans$create_project_structure_ans$delivery_dir, "norm_las")
      dir.create(to_dir, showWarnings = FALSE)
      # which files
      fls <- list.files(cloud2raster_ans$create_project_structure_ans$las_normalize_dir, full.names = T)
      # copy files
      q_file.copy <- purrr::quietly(file.copy)
      xxx <- fls %>%
        purrr::map(\(x) q_file.copy(from = x, to = to_dir))
      # delete temp files
      unlink(cloud2raster_ans$create_project_structure_ans$temp_dir, recursive = T)
    }

    ##### write foresttype raster data
    if( inherits(trees_type_rast, "SpatRaster") ){
      # write
      terra::writeRaster(
        trees_type_rast
        , paste0(cloud2raster_ans$create_project_structure_ans$delivery_dir, "/fia_foresttype_raster.tif")
        , overwrite = T
      )
    }

    ##### write biomass raster data
    if( inherits(trees_biomass_cell_landfire, "data.frame") ){
      # write
      write.csv(
        trees_biomass_cell_landfire
        , paste0(cloud2raster_ans$create_project_structure_ans$delivery_dir, "/stand_cell_data_landfire.csv")
        , row.names = F
      )
    }
    if( inherits(trees_biomass_cell_cruz, "data.frame") ){
      # write
      write.csv(
        trees_biomass_cell_cruz
        , paste0(cloud2raster_ans$create_project_structure_ans$delivery_dir, "/stand_cell_data_cruz.csv")
        , row.names = F
      )
    }

    # write settings and timer data
    xx11_fin <- Sys.time()
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
            , timer_trees_dbh_mins = difftime(xx6_trees_cbh, xx5_trees_dbh, units = c("mins")) %>%
              as.numeric()
            , timer_trees_cbh_mins = difftime(xx7_trees_type, xx6_trees_cbh, units = c("mins")) %>%
              as.numeric()
            , timer_trees_type_mins = difftime(xx8_trees_hmd, xx7_trees_type, units = c("mins")) %>%
              as.numeric()
            , timer_trees_hmd_mins = difftime(xx9_trees_biomass, xx8_trees_hmd, units = c("mins")) %>%
              as.numeric()
            , timer_trees_biomass_mins = difftime(xx10_return, xx9_trees_biomass, units = c("mins")) %>%
              as.numeric()
            , timer_write_data_mins = difftime(xx11_fin, xx10_return, units = c("mins")) %>%
              as.numeric()
            , timer_total_time_mins = difftime(xx11_fin, xx1_cloud2raster, units = c("mins")) %>%
              as.numeric()
            # settings
            , sttng_input_las_dir = input_las_dir[1]
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
            , sttng_dbh_model_regional = dbh_model_regional
            , sttng_dbh_model_local = dbh_model_local
            , sttng_estimate_dbh_from_cloud = estimate_dbh_from_cloud
            , sttng_estimate_tree_competition = estimate_tree_competition
            , sttng_competition_buffer_m = competition_buffer_m
            , sttng_competition_max_search_dist_m = competition_max_search_dist_m
            , sttng_estimate_tree_type = estimate_tree_type
            , sttng_type_max_search_dist_m = type_max_search_dist_m
            , sttng_estimate_tree_hmd = estimate_tree_hmd
            , sttng_hmd_tree_sample_n = hmd_tree_sample_n
            , sttng_hmd_tree_sample_prop = hmd_tree_sample_prop
            , sttng_hmd_estimate_missing_hmd = hmd_estimate_missing_hmd
            , sttng_estimate_biomass_method = which_biomass_methods %>% paste(collapse = ",")
            , sttng_biomass_max_crown_kg_per_m3 = biomass_max_crown_kg_per_m3
            , sttng_estimate_tree_cbh = estimate_tree_cbh
            , sttng_cbh_tree_sample_n = cbh_tree_sample_n
            , sttng_cbh_tree_sample_prop = cbh_tree_sample_prop
            , sttng_cbh_which_cbh = cbh_which_cbh
            , sttng_cbh_estimate_missing_cbh = cbh_estimate_missing_cbh
            , sttng_cbh_min_vhp_n = cbh_min_vhp_n
            , sttng_cbh_voxel_grain_size_m = cbh_voxel_grain_size_m
            , sttng_cbh_dist_btwn_bins_m = cbh_dist_btwn_bins_m
            , sttng_cbh_min_fuel_layer_ht_m = cbh_min_fuel_layer_ht_m
            , sttng_cbh_lad_pct_gap = cbh_lad_pct_gap
            , sttng_cbh_lad_pct_base = cbh_lad_pct_base
            , sttng_cbh_num_jump_steps = cbh_num_jump_steps
            , sttng_cbh_min_lad_pct = cbh_min_lad_pct
            , sttng_cbh_frst_layer_min_ht_m = cbh_frst_layer_min_ht_m
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
      , round(as.numeric(difftime(xx11_fin, xx1_cloud2raster, units = c("mins"))),2)
      , " minutes to process "
      , scales::comma(sum(cloud2raster_ans$chunk_las_catalog_ans$las_ctg@data$Number.of.point.records))
      , " points over an area of "
      , scales::comma(as.numeric(cloud2raster_ans$chunk_las_catalog_ans$las_ctg@data$geometry %>% sf::st_union() %>% sf::st_area())/10000,accuracy = 0.01)
      , " hectares"
    )
    #####################################
    # print errors
    #####################################
      if(!is.null(err_trees_competition)){
        message(paste0(
          "ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! in: trees_competition()"
          , "\n"
          , err_trees_competition
          , "\n..............try to run trees_competition() with updated parameter settings"
          , "\nusing `final_detected_tree_tops.gpkg` in the `point_cloud_processing_delivery` directory"
        ))
      }
      if(!is.null(err_treels_stem_dbh)){
        message(paste0(
          "ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! in: treels_stem_dbh()"
          , "\n"
          , err_treels_stem_dbh
          , "\n..............try to run treels_stem_dbh() with updated parameter settings"
          , "\nusing las files in `norm_las` in the `point_cloud_processing_delivery` directory"
        ))
      }
      if(!is.null(err_trees_dbh)){
        message(paste0(
          "ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! in: trees_dbh()"
          , "\n"
          , err_trees_dbh
          , "\n..............try to run trees_dbh() with updated parameter settings"
          , "\nusing `final_detected_tree_tops.gpkg` in the `point_cloud_processing_delivery` directory"
        ))
      }
      if(!is.null(err_trees_cbh)){
        message(paste0(
          "ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! in: trees_cbh()"
          , "\n"
          , err_trees_cbh
          , "\n..............try to run trees_cbh() with updated parameter settings"
          , "\nusing `final_detected_crowns.gpkg` and `norm_las` in the `point_cloud_processing_delivery` directory"
        ))
      }
      if(err_trees_cbh_sample==T){
        message(paste0(
          "WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! in: trees_cbh()"
          , "\n"
          , "The maximum number of trees to extract tree CBH using cloud2trees() is 20,000"
          , "\n..............try trees_cbh() if want to attempt to extract CBH for >20,000 trees"
          , "\nusing `final_detected_crowns.gpkg` and `norm_las` in the `point_cloud_processing_delivery` directory"
        ))
      }
      if(!is.null(err_trees_type)){
        message(paste0(
          "ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! in: trees_type()"
          , "\n"
          , err_trees_type
          , "\n..............try to run trees_type() with updated parameter settings"
          , "\nusing `final_detected_tree_tops.gpkg` in the `point_cloud_processing_delivery` directory"
        ))
      }
      if(!is.null(err_trees_hmd)){
        message(paste0(
          "ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! in: trees_hmd()"
          , "\n"
          , err_trees_hmd
          , "\n..............try to run trees_hmd() with updated parameter settings"
          , "\nusing `final_detected_crowns.gpkg` and `norm_las` in the `point_cloud_processing_delivery` directory"
        ))
      }
      if(err_trees_hmd_sample==T){
        message(paste0(
          "WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! in: trees_hmd()"
          , "\n"
          , "The maximum number of trees to extract tree HMD using cloud2trees() is 20,000"
          , "\n..............try trees_hmd() if want to attempt to extract HMD for >20,000 trees"
          , "\nusing `final_detected_crowns.gpkg` and `norm_las` in the `point_cloud_processing_delivery` directory"
        ))
      }
      if(!is.null(err_trees_biomass)){
        message(paste0(
          "ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! in: trees_biomass()"
          , "\n"
          , err_trees_biomass
          , "\n..............try to run trees_biomass() with updated parameter settings"
          , "\nusing `final_detected_tree_tops.gpkg` in the `point_cloud_processing_delivery` directory"
        ))
      }

    # return
    return(list(
      crowns_sf = crowns_sf_with_dbh
      , treetops_sf = treetops_sf_with_dbh
      , dtm_rast = cloud2raster_ans$dtm_rast
      , chm_rast = cloud2raster_ans$chm_rast
      , foresttype_rast = trees_type_rast
    ))
}
