#' @title Estimate CBH using tree crown polygons and normalized point cloud data
#'
#' @description
#' `trees_cbh()` uses the input tree crown polygons (e.g. as exported by [raster2trees()]) with the columns
#' `treeID` and `tree_height_m` to estimate tree CBH using height normalized point cloud data (e.g. as exported by [cloud2raster()]).
#'
#' CBH is extracted directly from the height normalized point cloud using the process outlined in Viedma et al. (2024) and implemented via [ladderfuelsr_cbh()].
#'
#' There are likely to be trees for which there is insufficient data in the point cloud to successfully estimate CBH. The user can elect to estimate missing CBH values which is accomplished via:
#'
#' * Attempt to extract CBH from the sample of trees elected by the user (`tree_sample_n`,`tree_sample_prop` parameter) using [ladderfuelsr_cbh()]
#' * Successfully extracted CBH trees become training data used to estimate the height-CBH allometry relationship that is spatially informed using the relative tree location compared to the training data
#' * The height and location predicting CBH model built from the point cloud training data is used to predict CBH for the non-training (i.e. missing CBH) data
#'
#' @param trees_poly must be one of the following that has required attributes `treeID` and `tree_height_m`:
#'   * `sf` class object with POLYGON geometry (see [sf::st_geometry_type()]). Recommended for smaller tree lists (e.g. <100k) that can fit in memory.
#'   * character vector with the path to a single or multiple spatial files that can be read by [sf::st_read()] and have POLYGON geometry. Recommended for large tree lists (e.g. 100k+) that might cause memory issues.
#'   * character with the path to a directory that has "final_detected_crowns*" files from [cloud2trees()] or [raster2trees()]. Recommended for large tree lists (e.g. 100k+) that might cause memory issues.
#'
#' @param norm_las character. a directory with nomalized las files, the path of a single .laz|.las file", -or- an object of class `LAS`.
#'   It is your responsibility to ensure that the point cloud is projected the same as the `trees_poly` data
#' @param tree_sample_n,tree_sample_prop numeric. Provide either `tree_sample_n`, the number of trees, or `tree_sample_prop`, the
#'   proportion of the trees to attempt to extract a CBH from the point cloud for.
#'   If neither are supplied, `tree_sample_n = 333` will be used. If both are supplied, `tree_sample_n` will be used.
#'   Increasing `tree_sample_prop` toward one (1) will increase the processing time, perhaps significantly depending on the number of trees in the `trees_poly` data.
#' @param which_cbh character. One of: "lowest"; "highest"; or "max_lad". See Viedma et al. (2024) reference.
#'   * "lowest" - Height of the CBH of the segmented tree based on the last distance found in its profile
#'   * "highest" - Height of the CBH of the segmented tree based on the maximum distance found in its profile
#'   * "max_lad" - Height of the CBH of the segmented tree based on the maximum LAD percentage
#' @param estimate_missing_cbh logical. even if the `tree_sample_prop` parameter is set to "1", it is not likely that CBH will be extracted successfully from every tree.
#'   Should the missing CBH values be estimated using the tree height and location information based on trees for which CBH is successfully extracted?
#' @param min_vhp_n numeric. the minimum number of vertical height profiles (VHPs) needed to estimate a CBH.
#' @param voxel_grain_size_m numeric. horizontal resolution (suggested 1 meter for lad profiles and 10 meters for LAI maps). See `grain.size` in [leafR::lad.voxels()]
#' @param dist_btwn_bins_m numeric. value for the actual height bin step (in meters). See `step` in [LadderFuelsR::get_gaps_fbhs()]
#' @param min_fuel_layer_ht_m numeric. value for the actual minimum base height (in meters). See `min_height` in [LadderFuelsR::get_gaps_fbhs()]
#' @param lad_pct_gap numeric. value of the percentile threshold used to identify gaps (default percentile 25th). See `perc_gap` in [LadderFuelsR::get_gaps_fbhs()]
#' @param lad_pct_base numeric. value of the percentile threshold used to identify fuels layers base height (default percentile 25th). See `perc_base` in [LadderFuelsR::get_gaps_fbhs()]
#' @param num_jump_steps numeric. value for the number of height bin steps that can be jumped to reshape fuels layers. See `number_steps` in [LadderFuelsR::get_real_fbh()]
#' @param min_lad_pct numeric. value for the minimum required LAD percentage in a fuel layer. See `threshold` in [LadderFuelsR::get_layers_lad()]
#' @param frst_layer_min_ht_m numeric. value for the depth height of the first fuel layer. If the first fuel layer has the maximum LAD and its depth is greater than the indicated value, then this fuel layer is considered as the CBH of the tree. On the contrary, if its depth is <= the value, the CBH with maximum LAD will be the second fuel layer, although it has not the maximum LAD. See `hdepth1_height` in [LadderFuelsR::get_cbh_metrics()]
#' @param force_same_crs logical. force the same crs between the point cloud and polygon if confident that data are in same projection.
#' data created by a `cloud2trees` pipeline (e.g. [cloud2raster()]) will always have the same projection even if not recognized by `lidR` functions
#' @param outfolder string. The path of a folder to write the model data to.
#'   Note, in the actual missing value estimation many RF models are estimated and model averaging is used.
#'   However, only the first estimated model is saved in this export which does not fully represent the process used to fill in missing values.
#'
#' @references
#' * [https://doi.org/10.1111/2041-210X.14427](https://doi.org/10.1111/2041-210X.14427)
#' Viedma, O., Silva, C. A., Moreno, J. M., & Hudak, A. T. (2024). LadderFuelsR: A new automated tool for vertical fuel continuity analysis and crown base height detection using light detection and ranging. Methods in Ecology and Evolution.
#' [https://github.com/olgaviedma/LadderFuelsR](https://github.com/olgaviedma/LadderFuelsR)
#'
#' * [https://doi.org/10.3390/rs11010092](https://doi.org/10.3390/rs11010092)
#' Almeida, D. R. A. D., Stark, S. C., Shao, G., Schietti, J., Nelson, B. W., Silva, C. A., ... & Brancalion, P. H. S. (2019). Optimizing the remote detection of tropical rainforest structure with airborne lidar: Leaf area profile sensitivity to pulse density and spatial sampling. Remote Sensing, 11(1), 92.
#' [https://github.com/DRAAlmeida/leafR](https://github.com/DRAAlmeida/leafR)
#'
#' @return Returns a spatial data frame of individual trees.
#'
#' @examples
#'  \dontrun{
#'  library(tidyverse)
#'  library(sf)
#'  # example tree crown polygons
#'  f <- system.file(package = "cloud2trees","extdata","crowns_poly.gpkg")
#'  crowns <- sf::st_read(f, quiet = T)
#'  # example normalized las files are in this directory
#'  norm_d <- system.file(package = "cloud2trees","extdata","norm_las")
#'  # now run the trees_cbh()
#'  trees_cbh_ans <- trees_cbh(
#'     trees_poly = crowns
#'     , norm_las = norm_d
#'     , tree_sample_n = 44
#'     , estimate_missing_cbh = T
#'     , force_same_crs = T
#'    )
#'  # what?
#'  trees_cbh_ans %>% class()
#'  trees_cbh_ans %>% dplyr::select(treeID,tidyselect::contains("cbh")) %>% dplyr::glimpse()
#'  # spatial polygons
#'  trees_cbh_ans %>% ggplot2::ggplot() +
#'    ggplot2::geom_sf(ggplot2::aes(fill=tree_cbh_m,color=is_training_cbh))
#'  # relationship between height and cbh
#'  trees_cbh_ans %>%
#'     ggplot2::ggplot(
#'       ggplot2::aes(x = tree_height_m, y = tree_cbh_m, color=is_training_cbh)
#'      ) +
#'     ggplot2::geom_point()
#'  # tabulate training data
#'  trees_cbh_ans %>%
#'    sf::st_drop_geometry() %>%
#'    dplyr::count(is_training_cbh)
#'  #### try a file list
#'  #### Recommended for large tree lists (e.g. 100k+) that might cause memory issues.
#'  # we'll split the crowns
#'  # as is done automatically for tree lists >250k by raster2trees() and cloud2trees()
#'  crowns <- crowns %>%
#'    dplyr::mutate(
#'      # makes 2 groups of data
#'      grp = ceiling(dplyr::row_number()/(dplyr::n()/2))
#'    )
#'  # make file names
#'  my_dir <- tempdir()
#'  fnm_1 <- file.path(my_dir, "crowns1.gpkg")
#'  fnm_2 <- file.path(my_dir, "crowns2.gpkg")
#'  fnm_1
#'  # write the data
#'  sf::st_write(crowns %>% dplyr::filter(grp==1), dsn = fnm_1, append = F) # grp 1
#'  sf::st_write(crowns %>% dplyr::filter(grp==2), dsn = fnm_2, append = F) # grp 2
#'  # try trees_cbh with our file list
#'  flist <- c(fnm_1,fnm_2)
#'  # now run the trees_cbh()
#'  trees_cbh_ans2 <- trees_cbh(
#'   trees_poly = flist
#'   , norm_las = norm_d
#'   , tree_sample_n = 44
#'   , estimate_missing_cbh = T
#'   , force_same_crs = T
#'  )
#'  # tabulate training data
#'  trees_cbh_ans2 %>%
#'    sf::st_drop_geometry() %>%
#'    dplyr::count(is_training_cbh)
#'  }
#' @export
#'
trees_cbh <- function(
  trees_poly
  , norm_las = NULL
  , tree_sample_n = NA
  , tree_sample_prop = NA
  , which_cbh = "lowest"
  , estimate_missing_cbh = TRUE
  , min_vhp_n = 3
  , voxel_grain_size_m = 1
  , dist_btwn_bins_m = 1
  , min_fuel_layer_ht_m = 1
  , lad_pct_gap = 25
  , lad_pct_base = 25
  , num_jump_steps = 1
  , min_lad_pct = 10
  , frst_layer_min_ht_m = 1
  , force_same_crs = F
  , outfolder = tempdir()
){
  # could move to parameters
  force_cbh_lte_ht <- T
  estimate_missing_max_n_training <- 25000
  #############################################
  # estimate cbh based on what's in trees_poly
  #############################################
  crowns_flist <- NULL # setting up for filling
  if(inherits(trees_poly, "sf")){
    #############################################
    # if trees_poly is sf
    #############################################
    cbh_df <- trees_cbh_sf(
      trees_poly = trees_poly
      , norm_las = norm_las
      , tree_sample_n = tree_sample_n
      , tree_sample_prop = tree_sample_prop
      , which_cbh = which_cbh
      , min_vhp_n = min_vhp_n
      , voxel_grain_size_m = voxel_grain_size_m
      , dist_btwn_bins_m = dist_btwn_bins_m
      , min_fuel_layer_ht_m = min_fuel_layer_ht_m
      , lad_pct_gap = lad_pct_gap
      , lad_pct_base = lad_pct_base
      , num_jump_steps = num_jump_steps
      , min_lad_pct = min_lad_pct
      , frst_layer_min_ht_m = frst_layer_min_ht_m
      , force_same_crs = force_same_crs
    )
  }else if(inherits(trees_poly, "character")){
    msg <- paste0(
      "If attempting to pass a list of files, the file list must:"
      , "\n   * be a vector of class character -AND-"
      , "\n   * be a directory that has final_detected_crowns* files from cloud2trees::cloud2trees() or cloud2trees::raster2trees()"
      , "\n   * -OR- be a vector of class character that includes spatial files that can be read by sf::st_read()"
    )
    trees_poly <- trees_poly %>% normalizePath() %>% unique()
    #############################################
    # if trees_poly is character
    #############################################
    # check if is dir and look for cloud2trees files in the dir
    if(
      length(trees_poly) == 1
      && dir.exists(trees_poly)
    ){
      #############################################
      # if trees_poly is a directory
      #############################################
      search_dir_final_detected_ans <- search_dir_final_detected(dir = trees_poly)
      crowns_flist <- search_dir_final_detected_ans$crowns_flist
      ttops_flist <- search_dir_final_detected_ans$ttops_flist
      if(is.null(crowns_flist)){
        stop(msg)
      }else if(length(crowns_flist)==1){ # if only one file then trees_cbh_sf() so that only read file once
        cbh_df <- trees_cbh_sf(
          trees_poly = crowns_flist
          , norm_las = norm_las
          , tree_sample_n = tree_sample_n
          , tree_sample_prop = tree_sample_prop
          , which_cbh = which_cbh
          , min_vhp_n = min_vhp_n
          , voxel_grain_size_m = voxel_grain_size_m
          , dist_btwn_bins_m = dist_btwn_bins_m
          , min_fuel_layer_ht_m = min_fuel_layer_ht_m
          , lad_pct_gap = lad_pct_gap
          , lad_pct_base = lad_pct_base
          , num_jump_steps = num_jump_steps
          , min_lad_pct = min_lad_pct
          , frst_layer_min_ht_m = frst_layer_min_ht_m
          , force_same_crs = force_same_crs
        )
      }else{
        cbh_df <- trees_cbh_flist(
          flist = trees_poly # just searches the directory again so that tree points are used for sample
          , norm_las = norm_las
          , tree_sample_n = tree_sample_n
          , tree_sample_prop = tree_sample_prop
          , which_cbh = which_cbh
          , min_vhp_n = min_vhp_n
          , voxel_grain_size_m = voxel_grain_size_m
          , dist_btwn_bins_m = dist_btwn_bins_m
          , min_fuel_layer_ht_m = min_fuel_layer_ht_m
          , lad_pct_gap = lad_pct_gap
          , lad_pct_base = lad_pct_base
          , num_jump_steps = num_jump_steps
          , min_lad_pct = min_lad_pct
          , frst_layer_min_ht_m = frst_layer_min_ht_m
          , force_same_crs = force_same_crs
        )
      }
    }else if(length(trees_poly)==1 && !file.exists(trees_poly)){
      #############################################
      # if trees_poly is a filename that doesn't exist
      #############################################
      stop(paste0(
        "could not find the file:"
        , "\n    "
        , trees_poly
      ))
    }else if(length(trees_poly)==1){
      #############################################
      # if trees_poly is a filename that does exist
      #############################################
      crowns_flist <- trees_poly
      cbh_df <- trees_cbh_sf(
        trees_poly = trees_poly
        , norm_las = norm_las
        , tree_sample_n = tree_sample_n
        , tree_sample_prop = tree_sample_prop
        , which_cbh = which_cbh
        , min_vhp_n = min_vhp_n
        , voxel_grain_size_m = voxel_grain_size_m
        , dist_btwn_bins_m = dist_btwn_bins_m
        , min_fuel_layer_ht_m = min_fuel_layer_ht_m
        , lad_pct_gap = lad_pct_gap
        , lad_pct_base = lad_pct_base
        , num_jump_steps = num_jump_steps
        , min_lad_pct = min_lad_pct
        , frst_layer_min_ht_m = frst_layer_min_ht_m
        , force_same_crs = force_same_crs
      )
    }else{
      #############################################
      # if trees_poly is a list of filenames
      #############################################
      crowns_flist <- trees_poly %>% normalizePath() %>% unique()
      cbh_df <- trees_cbh_flist(
        flist = crowns_flist
        , norm_las = norm_las
        , tree_sample_n = tree_sample_n
        , tree_sample_prop = tree_sample_prop
        , which_cbh = which_cbh
        , min_vhp_n = min_vhp_n
        , voxel_grain_size_m = voxel_grain_size_m
        , dist_btwn_bins_m = dist_btwn_bins_m
        , min_fuel_layer_ht_m = min_fuel_layer_ht_m
        , lad_pct_gap = lad_pct_gap
        , lad_pct_base = lad_pct_base
        , num_jump_steps = num_jump_steps
        , min_lad_pct = min_lad_pct
        , frst_layer_min_ht_m = frst_layer_min_ht_m
        , force_same_crs = force_same_crs
      )
    }
  }else{
    stop(paste0(
      "`trees_poly` data must be: "
      , "\n   * an object of class `sf` with only POLYGON type"
      , "\n   * -OR- a directory that has final_detected_crowns* files from cloud2trees::cloud2trees() or cloud2trees::raster2trees()"
      , "\n   * -OR- a vector of class character that includes spatial files that can be read by sf::st_read()"
    ))
  }

  #############################################
  # read cbh_df if it was a file list to get
  # trees where cbh extraction was successful
  #############################################
  if(
    inherits(cbh_df, "character")
    && (stringr::str_ends(cbh_df, ".*\\.csv$") %>% any())
  ){
    # read the output file(s)
    cbh_df <- stringr::str_subset(cbh_df, pattern = ".*\\.csv$") %>%
      readr::read_csv(progress = F, show_col_types = F) %>%
      dplyr::mutate(dplyr::across(
        .cols = c(tree_height_m, tree_cbh_m, tidyselect::ends_with("_zzz"))
        , .fns = as.numeric
      )) %>%
      dplyr::select(
        treeID, tree_height_m, tree_cbh_m, is_training_cbh, tidyselect::ends_with("_zzz")
      ) %>%
      dplyr::filter(is_training_cbh==T)
  }else if(inherits(cbh_df, "data.frame")){
    cbh_df <- cbh_df %>%
      dplyr::mutate(dplyr::across(
        .cols = c(tree_height_m, tree_cbh_m, tidyselect::ends_with("_zzz"))
        , .fns = as.numeric
      )) %>%
      dplyr::select(
        treeID, tree_height_m, tree_cbh_m, is_training_cbh, tidyselect::ends_with("_zzz")
      ) %>%
      dplyr::filter(is_training_cbh==T)
  }else{
    stop("error extracting CBH")
  }
  # ensure that there are enough data to estimate
  n_cbh <- nrow(cbh_df)

  #######################################################
  # build cbh model using training data
  # happens before we load trees_poly if inherits(crowns_flist, "character")
  # so that memory isn't already stuffed
  #######################################################
  n_cbh <- dplyr::coalesce(n_cbh,0)
  if(
    estimate_missing_cbh==T
    && n_cbh > 10
    && (names(cbh_df) %>% stringr::str_equal("tree_height_m") %>% any())
  ){
    # go for at least 50% of training data sampled
    # it's just that the individual model fits will be smaller
    ntimes_temp <- ((nrow(cbh_df)*0.5)/estimate_missing_max_n_training) %>%
      ceiling() %>%
      max(3) %>%
      min(50)
    # estimate
    cbh_mod <- rf_subsample_and_model_n_times(
      predictors = cbh_df %>% dplyr::select(-c(treeID,tree_cbh_m,is_training_cbh))
      , response = cbh_df$tree_cbh_m
      , mod_n_subsample = estimate_missing_max_n_training
      , mod_n_times = ntimes_temp
    )
    # write just the first model
    if(
      !is.null(cbh_mod)
      && length(cbh_mod)>0
    ){
      saveRDS(cbh_mod[[1]], file = file.path(normalizePath(outfolder), "cbh_height_model_estimates.rds"))
    }
  }else{
    cbh_mod <- NULL
  }

  #############################################
  # read trees_poly data to get full tree list
  #############################################
  if(inherits(trees_poly, "sf")){
    # get rid of columns we'll create
    trees_poly <- trees_poly %>%
      # throw in hey_xxxxxxxxxx to test it works if we include non-existant columns
      dplyr::select( -dplyr::any_of(c(
        "hey_xxxxxxxxxx"
        , "tree_cbh_m"
        , "is_training_cbh"
      )))
  }else if(
    !is.null(crowns_flist)
    && inherits(crowns_flist, "character")
  ){
    # if we've made it this far, the polygon data has already gone through the checks in trees_cbh_sf()
    # but we'll check again as it is quick
    # check_trees_poly will throw error if fails any checks
    check_trees_poly_ans <- crowns_flist %>% purrr::map(check_trees_poly)
    # read it to get the full list of tree polygons
    trees_poly <- crowns_flist %>%
      purrr::map(function(x){
        sf::st_read(
          dsn = x
          , quiet = T
        ) %>%
        # throw in hey_xxxxxxxxxx to test it works if we include non-existant columns
        dplyr::select( -dplyr::any_of(c(
          "hey_xxxxxxxxxx"
          , "tree_cbh_m"
          , "is_training_cbh"
        )))
      }) %>%
      dplyr::bind_rows()
  }else{
    stop("could not find tree crown polygon data")
  }

  #############################################
  # check for same class of treeID
  #############################################
  id_class <- class(trees_poly$treeID)[1]
  # cast treeID in original type
  if(!inherits(cbh_df$treeID, id_class)){
    if(id_class=="character"){
      cbh_df <- cbh_df %>%
        dplyr::mutate(treeID = as_character_safe(treeID))
    }
    if(id_class=="numeric"){
      cbh_df <- cbh_df %>%
        dplyr::mutate(treeID = as.numeric(treeID))
    }
  }

  #######################################################
  # fill missing cbh values
  #######################################################
  f <- trees_poly %>% names() %>% dplyr::coalesce("")
  n_cbh <- dplyr::coalesce(n_cbh,0)
  if(
    estimate_missing_cbh==T
    && !is.null(cbh_mod)
    && length(cbh_mod)>0
    && (names(trees_poly) %>% stringr::str_equal("tree_height_m") %>% any())
  ){
    # model the missing values
    predict_df <- trees_poly %>%
      dplyr::anti_join(cbh_df %>% dplyr::select(treeID), by = "treeID") %>%
      dplyr::select(treeID, tree_height_m) %>%
      make_spatial_predictors()
    # get predicted values
    predicted_temp <- rf_model_avg_predictions(mod_list = cbh_mod, predict_df = predict_df)
    # attach back to data
    predict_df <- predict_df %>%
      dplyr::mutate(predicted_zzz = predicted_temp$predicted)

    ## combine predicted data with training data for full data set
    trees_poly <- trees_poly %>%
      # join with training data
      dplyr::left_join(
        cbh_df %>%
          dplyr::select(treeID, tree_cbh_m, is_training_cbh)
        , by = "treeID"
      ) %>%
      dplyr::mutate(is_training_cbh = dplyr::coalesce(is_training_cbh, F)) %>%
      # join with predicted data estimates
      dplyr::left_join(
        predict_df %>%
          dplyr::select(treeID, predicted_zzz)
        , by = dplyr::join_by("treeID")
      ) %>%
      # clean up data
      dplyr::mutate(
        tree_cbh_m = dplyr::coalesce(tree_cbh_m, predicted_zzz)
      ) %>%
      dplyr::select(-predicted_zzz)

  }else if(n_cbh==0){
    message(paste0(
      "No CBH values extracted"
    ))
    return(
      trees_poly %>%
        dplyr::mutate(
          tree_cbh_m = as.numeric(NA)
          , is_training_cbh = as.logical(NA)
        )
    )
  }else if(estimate_missing_cbh==T){
    if(
      !(names(trees_poly) %>% stringr::str_equal("tree_height_m") %>% any())
    ){
      message(paste0(
        "`trees_poly` data must contain `tree_height_m` column to estimate CBH."
        , "\nSetting `estimate_missing_cbh=TRUE` requires this data."
        , "\nReturning CBH values extracted from cloud only."
      ))
    }else{
      message(paste0(
        "Insufficient data available to estimate missing CBH values."
        , "\nReturning CBH values extracted from cloud only."
      ))
    }

    ## combine predicted data with training data for full data set
    trees_poly <- trees_poly %>%
      # join with training data estimates
      dplyr::left_join(
        cbh_df %>%
          dplyr::filter(is_training_cbh==T) %>%
          dplyr::select(treeID, tree_cbh_m, is_training_cbh)
        , by = "treeID"
      ) %>%
      dplyr::mutate(is_training_cbh = dplyr::coalesce(is_training_cbh, F))
  }else{
    ## combine predicted data with training data for full data set
    trees_poly <- trees_poly %>%
      # join with training data estimates
      dplyr::left_join(
        cbh_df %>%
          dplyr::filter(is_training_cbh==T) %>%
          dplyr::select(treeID, tree_cbh_m, is_training_cbh)
        , by = "treeID"
      ) %>%
      dplyr::mutate(is_training_cbh = dplyr::coalesce(is_training_cbh, F))
  }

  ## prevent the CBH from being > the tree height
    if(
      force_cbh_lte_ht==T &&
      (names(trees_poly) %>% stringr::str_equal("tree_height_m") %>% any())
    ){
      # find the 95th percentile of height-cbh ratio
      max_ratio <- trees_poly %>%
        dplyr::filter(
          is_training_cbh==T
          & tree_cbh_m < tree_height_m
        ) %>%
        dplyr::mutate(ratio = tree_cbh_m/tree_height_m) %>%
        dplyr::pull(ratio) %>%
        stats::quantile(probs = 0.95)
      # update values
      trees_poly <- trees_poly %>%
        dplyr::mutate(
          # update training data where tree_cbh_m > tree_height_m
          is_training_cbh = dplyr::case_when(
            is_training_cbh==T & tree_cbh_m >= tree_height_m ~ FALSE
            , T ~ is_training_cbh
          )
          # update tree_cbh_m
          , tree_cbh_m = dplyr::case_when(
            is_training_cbh==F & tree_cbh_m/tree_height_m > max_ratio ~ max_ratio*tree_height_m
            , T ~ tree_cbh_m
          )
        )
    }
  # return
  return(trees_poly)

}
