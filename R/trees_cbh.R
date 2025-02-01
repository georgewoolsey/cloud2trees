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
#' @param trees_poly sf. A `sf` class object with POLYGON geometry (see [sf::st_geometry_type()]), the program will use the data "as-is" and only require the `treeID` and `tree_height_m` columns.
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
#'  # example tree crown polygons
#'  f <- system.file(package = "cloud2trees","extdata","crowns_poly.gpkg")
#'  crowns <- sf::st_read(f, quiet = T)
#'  # example normalized las files are in this directory
#'  norm_d <- system.file(package = "cloud2trees","extdata","norm_las")
#'  # now run the trees_cbh()
#'  trees_cbh_ans <- trees_cbh(
#'     trees_poly = crowns
#'     , norm_las = norm_d
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
){
  force_cbh_lte_ht <- T
  ##################################
  # check sample proportion
  ##################################
  if(
    is.na(as.numeric(tree_sample_n)) && is.na(as.numeric(tree_sample_prop))
  ){
    tree_sample_n <- 333
  }else if(
    !is.na(as.numeric(tree_sample_n)) && !is.na(as.numeric(tree_sample_prop))
  ){
    tree_sample_n <- dplyr::case_when(
      as.numeric(tree_sample_n)<=0 ~ 333
      , T ~ as.numeric(tree_sample_n)
    )
    tree_sample_prop <- NA
  }else if(
    is.na(as.numeric(tree_sample_n)) && !is.na(as.numeric(tree_sample_prop))
  ){
    tree_sample_prop <- dplyr::case_when(
      as.numeric(tree_sample_prop)<=0 ~ 0.5
      , as.numeric(tree_sample_prop)>1 ~ 1
      , T ~ as.numeric(tree_sample_prop)
    )
    tree_sample_n <- NA
  }else if(
    !is.na(as.numeric(tree_sample_n)) && is.na(as.numeric(tree_sample_prop))
  ){
    tree_sample_n <- dplyr::case_when(
      as.numeric(tree_sample_n)<=0 ~ 333
      , T ~ as.numeric(tree_sample_n)
    )
    tree_sample_prop <- NA
  }else{
    tree_sample_n <- 333
    tree_sample_prop <- NA
  }
  ##################################
  # check which cbh
  ##################################
    # clean it
    which_cbh <- dplyr::coalesce(which_cbh, "lowest")
    which_cbh <- ifelse(
      stringr::str_remove_all(which_cbh,"\\s") == ""
      , "lowest"
      , which_cbh
    )
    which_cbh <- tolower(which_cbh[1])
    # check it
    cbh_l <- c("max_lad", "highest", "lowest")
    if(!which_cbh %in% cbh_l){
      stop(paste0(
        "`which_cbh` must be one of:"
        , "\n"
        , paste(cbh_l, collapse = ", ")
      ))
    }
  ##################################
  # ensure that norm las data exists
  ##################################
  nlas_ctg <- check_las_data(norm_las)
  # set the lascatalog options
  if(inherits(nlas_ctg, "LAScatalog")){
    lidR::opt_progress(nlas_ctg) <- F
    lidR::opt_filter(nlas_ctg) <- "-drop_duplicates -drop_class 2 9 18" ## class 2 = ground; 9 = water; 18 = noise
    lidR::opt_select(nlas_ctg) <- "xyz0" # 0 enables all extra bytes to be loaded...possibly treeID
    lidR::opt_output_files(nlas_ctg) <- paste0(tempdir(), "/{*}_treed")
  }else if(inherits(nlas_ctg, "LAS")){
    stop(paste0(
      "`norm_las` should contain: a directory with nomalized las files,"
      ,"\n   the path of a single .laz|.las file,"
      , "\n   -or- an object of class `LAScatalog`"
    ))
    # nlas_ctg <- nlas_ctg %>%
    #   lidR::filter_poi(!Classification %in% c(2,9,18)) %>%
    #   lidR::filter_duplicates()
  }
  ##################################
  # ensure that treeID data exists
  ##################################
  f <- trees_poly %>% names() %>% dplyr::coalesce("")
  if(
    !(stringr::str_equal(f, "treeID") %>% any())
  ){
    stop(paste0(
      "`trees_poly` data must contain `treeID` column to estimate missing CBH values."
      , "\nProvide the `treeID` as a unique identifier of individual trees."
    ))
  }else{
    # check for duplicate treeID
    if(
      nrow(trees_poly) != length(unique(trees_poly$treeID))
    ){
      stop("Duplicates found in the treeID column. Please remove duplicates and try again.")
    }
  }
  # check for tree_height_m
  if(
    !(stringr::str_equal(f, "tree_height_m") %>% any())
  ){
    stop(paste0(
      "`trees_poly` data must contain `tree_height_m` column to estimate CBH."
      , "\nRename the height column if it exists and ensure it is in meters."
    ))
  }

  ##################################
  # ensure spatial polygon data
  ##################################
  sf_msg <- paste0(
      "`trees_poly` data must be an object of class `sf` with only POLYGON type."
      , "\nProvide an `sf` object and see `sf::st_geometry_type()`."
    )
  if(!inherits(trees_poly, "sf")){stop(sf_msg)}
  if( !(sf::st_is(trees_poly, type = c("POLYGON", "MULTIPOLYGON")) %>% all()) ){stop(sf_msg)}

  # get rid of columns we'll create
    trees_poly <- trees_poly %>%
      # throw in hey_xxxxxxxxxx to test it works if we include non-existant columns
      dplyr::select( -dplyr::any_of(c(
        "hey_xxxxxxxxxx"
        , "tree_cbh_m"
        , "is_training_cbh"
      )))

  ####################################################################
  # catalog apply
  ####################################################################
  # sample
    if(
      !is.na(tree_sample_prop)
      && tree_sample_prop<1
    ){
      samp_trees <- trees_poly %>%
        dplyr::slice_sample(
          prop = tree_sample_prop
        )
    }else if(
      !is.na(tree_sample_n)
      && tree_sample_n<nrow(trees_poly)
    ){
      samp_trees <- trees_poly %>%
        dplyr::slice_sample(
          n = tree_sample_n
        )
    }else{
      samp_trees <- trees_poly
    }
  ##################################
  # apply the ctg_leafr_for_ladderfuelsr function
  ##################################
  # simplify the polygons so that lidR::merge_spatial can be used
  simp_trees_poly <- simplify_multipolygon_crowns(samp_trees)
  # apply it
  output_temp <- lidR::catalog_apply(
    ctg = nlas_ctg
    , FUN = ctg_leafr_for_ladderfuelsr
    , .options = list(automerge = TRUE)
    # ctg_calc_tree_cbh options
    , poly_df = simp_trees_poly
    , force_crs = force_same_crs
    , voxel_grain_size_m = voxel_grain_size_m
  )

  ##################################
  # ladderfuelsr_cbh() to get cbh data by tree
  ##################################
  cbh_df <- output_ctg_to_ladderfuelsr_cbh(
    output = output_temp
    , min_vhp_n = min_vhp_n
    , voxel_grain_size_m = voxel_grain_size_m
    , dist_btwn_bins_m = dist_btwn_bins_m
    , min_fuel_layer_ht_m = min_fuel_layer_ht_m
    , lad_pct_gap = lad_pct_gap
    , lad_pct_base = lad_pct_base
    , num_jump_steps = num_jump_steps
    , min_lad_pct = min_lad_pct
    , frst_layer_min_ht_m = frst_layer_min_ht_m
  )
  ##################################
  # check the cbh data we got and use the cbh selected
  ##################################
  cbh_df <- clean_cbh_df(
    cbh_df = cbh_df
    , trees_poly = trees_poly
    , force_cbh_lte_ht = force_cbh_lte_ht
    , which_cbh = which_cbh
  )
  # ensure that there are enough data to estimate
  n_cbh <- nrow(cbh_df)

  ####################################################################
  # if zero CBH records were extracted
  ####################################################################
    # it is likely that the parameters entered caused an error in the `LadderFuelsR` workflow
    # the LadderFuelsR package commonly results in errors if the parameters don't align
    # i haven't been able to find a pattern other than:
    # * ERROR if : min_fuel_layer_ht_m == dist_btwn_bins_m (except 0.5,1.5,2.5)??
    # * ERROR if : dist_btwn_bins_m - min_fuel_layer_ht_m >= 1.5
    # ............. try to run it with settings that are close and maybe work better
    if(dplyr::coalesce(n_cbh,0)==0){
      ########################
      # new_min_fuel_layer_ht_m
      new_min_fuel_layer_ht_m <- dplyr::case_when(
        dist_btwn_bins_m - min_fuel_layer_ht_m >= 1.5 ~ min_fuel_layer_ht_m
        , min_fuel_layer_ht_m == dist_btwn_bins_m ~ min_fuel_layer_ht_m
        , T ~ 1 # default that has worked with every try
      )
      # new_dist_btwn_bins_m
      new_dist_btwn_bins_m <- dplyr::case_when(
        dist_btwn_bins_m - min_fuel_layer_ht_m >= 1.5 ~ min_fuel_layer_ht_m + 1
        , (min_fuel_layer_ht_m == dist_btwn_bins_m) & min_fuel_layer_ht_m >= 1 ~ min_fuel_layer_ht_m - 0.5
        , min_fuel_layer_ht_m == dist_btwn_bins_m ~ min_fuel_layer_ht_m + 0.5
        , T ~ 0.5 # default that has worked with every try
      )
      # default some others
      new_min_vhp_n <- 3 # default
      new_voxel_grain_size_m <- 2 # default
      new_min_lad_pct <- 10 # default
      # message
      message(paste0(
        "No CBH values extracted from point cloud with supplied parameters."
        , "\nAttempting to update parameters..."
        , "\nmin_fuel_layer_ht_m from: ", min_fuel_layer_ht_m, " --> to: ", new_min_fuel_layer_ht_m
        , "\ndist_btwn_bins_m from: ", dist_btwn_bins_m, " --> to: ", new_dist_btwn_bins_m
        , "\nmin_vhp_n from: ", min_vhp_n, " --> to: ", new_min_vhp_n
        , "\nvoxel_grain_size_m from: ", voxel_grain_size_m, " --> to: ", new_voxel_grain_size_m
        , "\nmin_lad_pct from: ", min_lad_pct, " --> to: ", new_min_lad_pct
      ))
      ##################################
      # apply the ctg_leafr_for_ladderfuelsr function
      ##################################
      # apply it
      if(new_voxel_grain_size_m!=voxel_grain_size_m){
        output_temp <- lidR::catalog_apply(
          ctg = nlas_ctg
          , FUN = ctg_leafr_for_ladderfuelsr
          , .options = list(automerge = TRUE)
          # ctg_calc_tree_cbh options
          , poly_df = simp_trees_poly
          , force_crs = force_same_crs
          , voxel_grain_size_m = new_voxel_grain_size_m
        )
      }

      ##################################
      # ladderfuelsr_cbh() to get cbh data by tree
      ##################################
      cbh_df <- output_ctg_to_ladderfuelsr_cbh(
        output = output_temp
        , min_vhp_n = new_min_vhp_n
        , voxel_grain_size_m = new_voxel_grain_size_m
        , dist_btwn_bins_m = new_dist_btwn_bins_m
        , min_fuel_layer_ht_m = new_min_fuel_layer_ht_m
        , lad_pct_gap = lad_pct_gap
        , lad_pct_base = lad_pct_base
        , num_jump_steps = num_jump_steps
        , min_lad_pct = new_min_lad_pct
        , frst_layer_min_ht_m = frst_layer_min_ht_m
      )
      ##################################
      # check the cbh data we got and use the cbh selected
      ##################################
      cbh_df <- clean_cbh_df(
        cbh_df = cbh_df
        , trees_poly = trees_poly
        , force_cbh_lte_ht = force_cbh_lte_ht
        , which_cbh = which_cbh
      )
      # ensure that there are enough data to estimate
      n_cbh <- nrow(cbh_df)

    }

  #############################################
  # check for estimate missing
  #############################################
  # ensure that tree height data exists
  f <- trees_poly %>% names() %>% dplyr::coalesce("")
  if(
    estimate_missing_cbh==T
    && n_cbh > 10
    && (names(trees_poly) %>% stringr::str_equal("tree_height_m") %>% any())
  ){
    # add x,y to data
    mod_df <- trees_poly %>%
      dplyr::left_join(
        cbh_df %>%
          dplyr::select(treeID, tree_cbh_m, is_training_cbh)
        , by = "treeID"
      ) %>%
      dplyr::mutate(is_training_cbh = dplyr::coalesce(is_training_cbh, F)) %>%
      dplyr::select(treeID, tree_height_m, tree_cbh_m, is_training_cbh) %>%
      dplyr::mutate(crown_area_zzz = sf::st_area(.) %>% as.numeric()) %>%
      sf::st_centroid() %>%
      dplyr::mutate(
        tree_xxx = sf::st_coordinates(.)[,1]
        , tree_yyy = sf::st_coordinates(.)[,2]
        , tree_height_m = as.numeric(tree_height_m)
        , tree_cbh_m = as.numeric(tree_cbh_m)
      ) %>%
      sf::st_drop_geometry()
    # training versus predict data
    training_df <- mod_df %>% dplyr::filter(is_training_cbh==T) %>% dplyr::select(-is_training_cbh)
    predict_df <- mod_df %>% dplyr::filter(is_training_cbh==F) %>% dplyr::select(-is_training_cbh)

    ### tuning RF model
      # If we are interested with just starting out and tuning the mtry parameter
      # we can use randomForest::tuneRF for a quick and easy tuning assessment.
      # tuneRf will start at a value of mtry that you supply and increase by a
      # certain step factor until the OOB error stops improving be a specified amount.
      # quiet this
      quiet_tuneRF <- purrr::quietly(randomForest::tuneRF)
      # run it
      rf_tune_temp <- quiet_tuneRF(
        # randomForest::tuneRF(
        y = training_df$tree_cbh_m
        , x = training_df %>% dplyr::select(-c(treeID,tree_cbh_m))
        , stepFactor = 0.5
        , ntreeTry = 200
        , mtryStart = 0.5
        , improve = 0.01
        , plot = F
        , trace = F
      )
      # just get the result
      rf_tune_temp <- rf_tune_temp$result

      # Extract the optimal mtry value
      optimal_mtry <- rf_tune_temp %>%
        dplyr::as_tibble() %>%
        dplyr::filter(OOBError==min(OOBError)) %>%
        dplyr::filter(dplyr::row_number() == 1) %>%
        dplyr::pull(mtry)
      # ensure that the mtry value is not greater than the number of predictors
      optimal_mtry <- min(
        optimal_mtry
        , ncol(
          training_df %>% dplyr::select(-c(treeID,tree_cbh_m))
        )
      )

      ### Run a randomForest model to predict CBH using various crown predictors
      # quiet this
      quiet_rf <- purrr::quietly(randomForest::randomForest)
      # run it
      cbh_mod <- quiet_rf(
        y = training_df$tree_cbh_m
        , x = training_df %>% dplyr::select(-c(treeID,tree_cbh_m))
        , mtry = optimal_mtry
        , na.action = na.omit
      )
      # just get the result
      cbh_mod <- cbh_mod$result

    # # model
    # cbh_mod <- stats::lm(
    #   formula = tree_cbh_m ~ tree_xxx + tree_yyy + tree_xxx:tree_yyy + tree_height_m + crown_area_zzz
    #   , data = training_df
    # )

    # predict missing
    predicted_cbh_temp <- predict(
        cbh_mod
        , predict_df %>% dplyr::select(-c(treeID,tree_cbh_m))
      ) %>%
      dplyr::as_tibble() %>%
      dplyr::pull(1)

    ## combine predicted data with training data for full data set
    trees_poly <- trees_poly %>%
      # join with training data estimates
      dplyr::left_join(
        cbh_df %>%
          dplyr::filter(is_training_cbh==T) %>%
          dplyr::select(treeID, tree_cbh_m, is_training_cbh)
        , by = "treeID"
      ) %>%
      dplyr::mutate(is_training_cbh = dplyr::coalesce(is_training_cbh, F)) %>%
      # join with predicted data estimates
      dplyr::left_join(
        predict_df %>%
          dplyr::mutate(
            predicted_cbh = predicted_cbh_temp
          ) %>%
          dplyr::select(treeID, predicted_cbh)
        , by = dplyr::join_by("treeID")
      ) %>%
      # clean up data
      dplyr::mutate(
        tree_cbh_m = dplyr::coalesce(tree_cbh_m, predicted_cbh)
      ) %>%
      dplyr::select(-predicted_cbh)

    ## prevent the CBH from being > the tree height
      # find the 95th percentile of height-cbh ratio
      max_ratio <- cbh_df %>%
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
  # return
  return(trees_poly)

}

#####################################################
#####################################################
# intermediate functions
#####################################################
#####################################################
####################################
## function to clip the point cloud to a polygon
## and run it through:
## leafr_for_ladderfuelsr() which is the only
## step in the process that uses the las
## and returns a data.frame
## to pass to ladderfuelsr_cbh()
####################################
ctg_leafr_for_ladderfuelsr <- function(
  chunk
  , poly_df
  , force_crs = F
  , voxel_grain_size_m = 1
){
  las <- lidR::readLAS(chunk)
  if(lidR::is.empty(las)){return(NULL)}
  # check for treeID in las already
  if(names(las@data) %>% stringr::str_equal("treeID") %>% any()){
    nlas_tree <- las
  }else{
    # attach treeID
    nlas_tree <- polygon_attribute_to_las(
      las
      , simplify_multipolygon_crowns(poly_df)
      , attribute = "treeID"
      , force_crs = force_crs
    )
  }
  # get the lad profile for each treeID
  safe_leafr_for_ladderfuelsr <- purrr::safely(leafr_for_ladderfuelsr)
  lad_profile <- safe_leafr_for_ladderfuelsr(
      nlas_tree
      , voxel_grain_size_m = voxel_grain_size_m
      , k = 1
      , attribute = "treeID"
      , relative = F
    )
  # just get the result
  lad_profile <- lad_profile$result
  # return
  return(lad_profile)
}
####################################
## function to handle the lidR::catalog_apply
## output and run it through:
## ladderfuelsr_cbh() and pull out cbh values
####################################
output_ctg_to_ladderfuelsr_cbh <- function(
  output = NULL
  , min_vhp_n = 3
  , voxel_grain_size_m = 1
  , dist_btwn_bins_m = 1
  , min_fuel_layer_ht_m = 1
  , lad_pct_gap = 25
  , lad_pct_base = 25
  , num_jump_steps = 1
  , min_lad_pct = 10
  , frst_layer_min_ht_m = 1
) {
  # make sure output is a readable file
  if(
    inherits(output, "character")
    && (stringr::str_ends(output, ".*\\.(txt|csv)$") %>% any())
  ){
    # read the output file(s)
    lad_profile <- stringr::str_subset(output, pattern = ".*\\.(txt|csv)$") %>%
      purrr::map(\(x) readr::read_delim(
        file = x, progress = F, show_col_types = F
      )) %>%
      dplyr::bind_rows() %>%
      # for trees on many tiles keep row with most points
      dplyr::group_by(treeID) %>%
      dplyr::filter(total_pulses == max(total_pulses)) %>%
      # if one tree has multiple cases with the same total pulses
      # need to get row unique by treeID, height
      # just take the first record
      # which is the first tile processed with the tree
      dplyr::group_by(treeID, height) %>%
      dplyr::summarise(dplyr::across(
         .cols = dplyr::everything()
         , .fns = dplyr::first
      )) %>%
      dplyr::ungroup()
  }else{
    return(NULL)
  }
  if(nrow(lad_profile)<1){return(NULL)}

  # force treeID to numeric
  if(!inherits(lad_profile$treeID, "numeric")){
    # get numeric
    lad_profile <- lad_profile %>%
      # make numeric
      dplyr::inner_join(
        lad_profile %>%
          dplyr::distinct(treeID) %>%
          dplyr::mutate(id = dplyr::row_number())
        , by = "treeID"
      ) %>%
      dplyr::mutate(
        treeID_bu=treeID
        , treeID = id
      ) %>%
      dplyr::select(-c(id)) %>%
      dplyr::relocate(treeID)
  }
  # dplyr::glimpse(lad_profile)
  # extract the CBH using ladderfuelsr_cbh()
  # we can map over multiple trees
  # quiet this function ... which should not issue an error
  # since all LadderFuelsR were handled with purrr::safely
  quiet_ladderfuelsr_cbh <- purrr::quietly(ladderfuelsr_cbh)
  cbh_df <- lad_profile$treeID %>%
    unique() %>%
    purrr::map(\(x)
        quiet_ladderfuelsr_cbh(
        # ladderfuelsr_cbh(
          lad_profile_df = lad_profile
          , treeID = x
          , min_vhp_n = min_vhp_n
          , voxel_grain_size_m = voxel_grain_size_m
          , dist_btwn_bins_m = dist_btwn_bins_m
          , min_fuel_layer_ht_m = min_fuel_layer_ht_m
          , lad_pct_gap = lad_pct_gap
          , lad_pct_base = lad_pct_base
          , num_jump_steps = num_jump_steps
          , min_lad_pct = min_lad_pct
          , frst_layer_min_ht_m = frst_layer_min_ht_m
        ) %>%
        purrr::pluck("result") %>% ## b/c purrr::quietly
        purrr::pluck("cbh_metrics")
      , .progress = "extracting CBH"
    ) %>%
    dplyr::bind_rows()
  # dplyr::glimpse(cbh_df)

  if(nrow(cbh_df)<1){return(NULL)}

  # clean cbh data
  cbh_df <- cbh_df %>%
    dplyr::mutate(
      cbh_maxlad_height_m = maxlad_Hcbh
      , cbh_max_height_m = max_Hcbh
      , cbh_last_height_m = last_Hcbh
    ) %>%
    dplyr::select(treeID, tidyselect::starts_with("cbh_")) %>%
    # just make sure that we didn't get multiple cbh by tree back
    dplyr::group_by(treeID) %>%
    dplyr::summarise(dplyr::across(
      tidyselect::starts_with("cbh_")
      , ~ mean(., na.rm = T)
    )) %>%
    dplyr::ungroup()

  # add the number of points (pulses) in the point cloud
  # and update treeID if made backup
  if(names(lad_profile) %>% stringr::str_equal("treeID_bu") %>% any()){
    cbh_df <- cbh_df %>%
      dplyr::mutate(treeID=as.character(treeID)) %>%
      dplyr::inner_join(
        lad_profile %>% dplyr::distinct(treeID_bu, treeID, total_pulses) %>%
          dplyr::mutate(treeID=as.character(treeID))
        , by = "treeID"
      ) %>%
      dplyr::mutate(
        treeID = treeID_bu
      ) %>%
      dplyr::select(-treeID_bu) %>%
      dplyr::relocate(treeID)
  }else{ # treeID was already numeric
    cbh_df <- cbh_df %>%
      dplyr::mutate(treeID=as.numeric(treeID)) %>% # LadderfuelsR turns treeID to factor
      dplyr::inner_join(
        lad_profile %>%
          dplyr::distinct(treeID, total_pulses) %>%
          dplyr::mutate(treeID=as.numeric(treeID))
        , by = "treeID"
      ) %>%
      dplyr::relocate(treeID)
  }

  # dplyr::glimpse(cbh_df)
  # return
  return(cbh_df)
}

# cbh_df <- output_ctg_to_ladderfuelsr_cbh(output)
# cbh_df %>% dplyr::glimpse()
# identical(
#   cbh_df %>% dplyr::distinct(treeID) %>% nrow
#   , cbh_df %>% nrow()
# )
# cbh_df %>%
#   tidyr::pivot_longer(
#     cols = tidyselect::starts_with("cbh_")
#   ) %>%
#   ggplot2::ggplot(
#     mapping = ggplot2::aes(x=value, color=name, fill=name)
#   ) +
#   ggplot2::geom_density() +
#   ggplot2::facet_grid(rows = dplyr::vars(name)) +
#   ggplot2::theme_light()

####################################
## function to check the
## return data from
## output_ctg_to_ladderfuelsr_cbh
####################################
clean_cbh_df <- function(cbh_df = NULL, trees_poly, force_cbh_lte_ht, which_cbh) {
  if(
    inherits(cbh_df, "data.frame")
    && dplyr::coalesce(nrow(cbh_df),0)>0
  ){
    # read the output file(s)
    cbh_df <- cbh_df %>%
      # for trees on many tiles keep row with most points
      dplyr::group_by(treeID) %>%
      dplyr::filter(total_pulses == max(total_pulses)) %>%
      dplyr::summarise(
        # if still duplicates
        cbh_maxlad_height_m = max(cbh_maxlad_height_m, na.rm = T)
        , cbh_max_height_m = max(cbh_max_height_m, na.rm = T)
        , cbh_last_height_m = min(cbh_last_height_m, na.rm = T)
      ) %>%
      dplyr::ungroup()
    # pick a cbh
    if(which_cbh == "max_lad"){
      cbh_df <- cbh_df %>%
        dplyr::mutate(
          tree_cbh_m = cbh_maxlad_height_m
        )
    }else if(which_cbh == "highest"){
      cbh_df <- cbh_df %>%
        dplyr::mutate(
          tree_cbh_m = cbh_max_height_m
        )
    }else{
      cbh_df <- cbh_df %>%
        dplyr::mutate(
          tree_cbh_m = cbh_last_height_m
        )
    }
  }else{
    # blank
    cbh_df <- dplyr::tibble(NULL)
  }

  # check force_cbh_lte_ht
  if(
    force_cbh_lte_ht==T &&
    (names(trees_poly) %>% stringr::str_equal("tree_height_m") %>% any()) &&
    nrow(cbh_df)>0
  ){
    # filter based on cbh vs ht
    cbh_df <- cbh_df %>%
      dplyr::inner_join(
        trees_poly %>%
          sf::st_drop_geometry() %>%
          dplyr::select(treeID,tree_height_m)
        , by = "treeID"
      ) %>%
      dplyr::filter(
        !is.na(tree_cbh_m)
        & tree_cbh_m < tree_height_m
      ) %>%
      dplyr::mutate(is_training_cbh=T)
  }else if(nrow(cbh_df)>0){
    # filter cbh
    cbh_df <- cbh_df %>%
      dplyr::filter(
        !is.na(tree_cbh_m)
      ) %>%
      dplyr::mutate(is_training_cbh=T)
  }

  # return
  return(cbh_df)
}
