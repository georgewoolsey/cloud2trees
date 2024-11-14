#' @title Estimate DBH for a tree list based on height
#'
#' @description
#' `trees_cbh()` uses the input tree crown polygons (e.g. as exported by [raster2trees()]) with the columns
#' `treeID` and `tree_height_m` to estimate tree CBH based on a normalized point cloud (e.g. as exported by [raster2trees()]).
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
#'   If neither are supplied, `tree_sample_n = 500` will be used. If both are supplied, `tree_sample_n` will be used.
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
#'  # example tree list
#'  }
#' @export
#'
trees_cbh <- function(
  trees_poly
  , norm_las = NA
  , tree_sample_n = NA
  , tree_sample_prop = NA
  , which_cbh = "lowest"
  , estimate_missing_cbh = FALSE
  , min_vhp_n = 4
  , voxel_grain_size_m = 2
  , dist_btwn_bins_m = 1
  , min_fuel_layer_ht_m = 1
  , lad_pct_gap = 25
  , lad_pct_base = 25
  , num_jump_steps = 1
  , min_lad_pct = 10
  , frst_layer_min_ht_m = 1
){
  ##################################
  # check sample proportion
  ##################################
  if(
    is.na(as.numeric(tree_sample_n)) && is.na(as.numeric(tree_sample_prop))
  ){
    tree_sample_n <- 500
  }else if(
    !is.na(as.numeric(tree_sample_n)) && !is.na(as.numeric(tree_sample_prop))
  ){
    tree_sample_n <- dplyr::case_when(
      as.numeric(tree_sample_n)<=0 ~ 500
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
      as.numeric(tree_sample_n)<=0 ~ 500
      , T ~ as.numeric(tree_sample_n)
    )
    tree_sample_prop <- NA
  }else{
    tree_sample_n <- 500
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
  nlas_msg <- paste0(
    "`norm_las` must contain a directory with nomalized las files, the path of a .laz|.las file"
    , "\n, -or- an object of class `LAS`. Please update the `norm_las` parameter."
  )
  if(is.na(norm_las)){stop(nlas_msg)}
  if(inherits(norm_las, "character")){
    if(!stringr::str_ends(norm_las, ".*\\.(laz|las)$")){
      # try to read directory for las files
      fls <- list.files(normalizePath(norm_las), pattern = ".*\\.(laz|las)$", full.names = TRUE)
      # stop it if no files
      if(length(fls)<1){stop(nlas_msg)}
      # read it
      nlas_ctg <- lidR::readLAScatalog(fls)
      # turn of lidR progress
      lidR::opt_progress(nlas_ctg) <- F
    }else if(stringr::str_ends(norm_las, ".*\\.(laz|las)$")){
      # read it
      nlas_ctg <- lidR::readLAScatalog(norm_las)
      # turn of lidR progress
      lidR::opt_progress(nlas_ctg) <- F
    }else{
      stop(nlas_msg)
    }
  }else if(inherits(norm_las, "LAS")){
    nlas_ctg <- norm_las
  }else{
    stop(nlas_msg)
  }
  ##################################
  # ensure that treeID data exists
  ##################################
  f <- trees_poly %>% names()
  if(length(f)==0){f <- ""}
  if(
    max(grepl("treeID", f))==0
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
    # ensure that treeID is numeric
    # generate a treeID index because it needs to be numeric for LadderFuelsR
    trees_poly <- trees_poly %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        treeID_backup = treeID
        , treeID = dplyr::row_number()
      ) %>%
      dplyr::relocate(treeID)
  }

  ##################################
  # ensure spatial polygon data
  ##################################
  sf_msg <- paste0(
      "`trees_poly` data must be an object of class `sf` with only POLYGON type."
      , "\nProvide an `sf` object and see `sf::st_geometry_type()`."
    )
  if(!inherits(trees_poly, "sf")){stop(sf_msg)}
  if( min(sf::st_is(trees_poly, type = c("POLYGON", "MULTIPOLYGON"))) == 0 ){stop(sf_msg)}

  ##################################
  # ensure the las and sf are same projection
  ##################################
  # get crs
    crs_las <- sf::st_crs(nlas_ctg)
    crs_poly <- sf::st_crs(trees_poly)
  # test equal epsg
    if(
      is.na(crs_las$epsg) |
      is.na(crs_poly$epsg) |
      crs_las$epsg != crs_poly$epsg
    ){
      # try to pull the epsg another way
      # get_horizontal_crs is defined in chunk_las_catalog.R
      n_crs <- get_horizontal_crs(nlas_ctg)
      if(
        is.na(n_crs) |
        n_crs$epsg != crs_poly$epsg
      ){
        stop("The `trees_poly` and `norm_las` data have differing CRS projections. Please see `sf::st_crs()` and ensure compatibility.")
      }
    }

  # get rid of columns we'll create
    trees_poly <- trees_poly %>%
      # throw in hey_xxxxxxxxxx to test it works if we include non-existant columns
      dplyr::select( -dplyr::any_of(c(
        "hey_xxxxxxxxxx"
        , "tree_cbh_m"
        , "is_training_cbh"
      )))

  ####################################################################
  # map over ladderfuelsr_cbh function
  ####################################################################
  # is this ladderfulesr thing even possible without errors even after safe running it in call_ladderfuelsr_cbh?????
  poss_call_ladderfuelsr_cbh <- purrr::possibly(call_ladderfuelsr_cbh, otherwise = NULL, quiet = T)
  # sample
    if(!is.na(tree_sample_prop)){
      samp_trees <- trees_poly %>%
        dplyr::slice_sample(
          prop = tree_sample_prop
        ) %>%
        dplyr::pull(treeID)
    }else{
      samp_trees <- trees_poly %>%
        dplyr::slice_sample(
          n = tree_sample_n
        ) %>%
        dplyr::pull(treeID)
    }
  # try it
  cbh_df <- samp_trees %>%
    purrr::map(\(x) poss_call_ladderfuelsr_cbh(
        id = x
        , poly_df = trees_poly
        , nlas = nlas_ctg
        , my_min_vhp_n = min_vhp_n
        , my_voxel_grain_size_m = voxel_grain_size_m
        , my_dist_btwn_bins_m = dist_btwn_bins_m
        , my_min_fuel_layer_ht_m = min_fuel_layer_ht_m
        , my_lad_pct_gap = lad_pct_gap
        , my_lad_pct_base = lad_pct_base
        , my_num_jump_steps = num_jump_steps
        , my_min_lad_pct = min_lad_pct
        , my_frst_layer_min_ht_m = frst_layer_min_ht_m
      )
      , .progress = "extracting CBH"
    ) %>%
    dplyr::bind_rows() %>%
    # get rid of treeID
    dplyr::mutate(
      treeID = treeID_backup
    ) %>%
    dplyr::select(-treeID_backup) %>%
    dplyr::relocate(treeID) %>%
    sf::st_drop_geometry()

  # pick a cbh
  if(which_cbh == "max_lad"){
    cbh_df <- cbh_df %>%
      dplyr::mutate(
        tree_cbh_m = cbh_maxlad_height_m
        , is_training_cbh = !is.na(cbh_maxlad_height_m)
      )
  }else if(which_cbh == "highest"){
    cbh_df <- cbh_df %>%
      dplyr::mutate(
        tree_cbh_m = cbh_max_height_m
        , is_training_cbh = !is.na(cbh_max_height_m)
      )
  }else{
    cbh_df <- cbh_df %>%
      dplyr::mutate(
        tree_cbh_m = cbh_last_height_m
        , is_training_cbh = !is.na(cbh_last_height_m)
      )
  }

  # ensure that there are enough data to estimate
  n_cbh <- cbh_df %>% dplyr::filter(is_training_cbh==T) %>% nrow()

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
      # new_min_vhp_n
      new_min_vhp_n <- 4 # default
      # message
      message(paste0(
        "No CBH values extracted from point cloud with supplied parameters."
        , "\nAttempting to update parameters..."
        , "\nmin_fuel_layer_ht_m from: ", min_fuel_layer_ht_m, " --> to: ", new_min_fuel_layer_ht_m
        , "\ndist_btwn_bins_m from: ", dist_btwn_bins_m, " --> to: ", new_dist_btwn_bins_m
        , "\nmin_vhp_n from: ", min_vhp_n, " --> to: ", new_min_vhp_n
      ))
      ####################################################################
      # map over ladderfuelsr_cbh function
      ####################################################################
        cbh_df <- samp_trees %>%
          purrr::map(\(x) poss_call_ladderfuelsr_cbh(
              id = x
              , poly_df = trees_poly
              , nlas = nlas_ctg
              , my_min_vhp_n = new_min_vhp_n
              , my_voxel_grain_size_m = voxel_grain_size_m
              , my_dist_btwn_bins_m = new_dist_btwn_bins_m
              , my_min_fuel_layer_ht_m = new_min_fuel_layer_ht_m
              , my_lad_pct_gap = lad_pct_gap
              , my_lad_pct_base = lad_pct_base
              , my_num_jump_steps = num_jump_steps
              , my_min_lad_pct = min_lad_pct
              , my_frst_layer_min_ht_m = frst_layer_min_ht_m
            )
            , .progress = "extracting CBH"
          ) %>%
          dplyr::bind_rows() %>%
          # get rid of treeID
          dplyr::mutate(
            treeID = treeID_backup
          ) %>%
          dplyr::select(-treeID_backup) %>%
          dplyr::relocate(treeID) %>%
          sf::st_drop_geometry()

        # pick a cbh
        if(which_cbh == "max_lad"){
          cbh_df <- cbh_df %>%
            dplyr::mutate(
              tree_cbh_m = cbh_maxlad_height_m
              , is_training_cbh = !is.na(cbh_maxlad_height_m)
            )
        }else if(which_cbh == "highest"){
          cbh_df <- cbh_df %>%
            dplyr::mutate(
              tree_cbh_m = cbh_max_height_m
              , is_training_cbh = !is.na(cbh_max_height_m)
            )
        }else{
          cbh_df <- cbh_df %>%
            dplyr::mutate(
              tree_cbh_m = cbh_last_height_m
              , is_training_cbh = !is.na(cbh_last_height_m)
            )
        }

        # ensure that there are enough data to estimate
        n_cbh <- cbh_df %>% dplyr::filter(is_training_cbh==T) %>% nrow()
    }

  #############################################
  # check for estimate missing
  #############################################
  # ensure that tree height data exists
  f <- trees_poly %>% names()
  if(length(f)==0){f <- ""}

  # get rid of treeID backup
  trees_poly <- trees_poly %>%
    dplyr::mutate(treeID = treeID_backup) %>%
    dplyr::select(-treeID_backup) %>%
    dplyr::relocate(treeID)

  if(
    estimate_missing_cbh==T
    & n_cbh > 10
    & max(grepl("tree_height_m", f))==1
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
      sf::st_centroid() %>%
      dplyr::mutate(
        tree_xxx = sf::st_coordinates(.)[,1]
        , tree_yyy = sf::st_coordinates(.)[,2]
        , crown_area_zzz = sf::st_area(.)
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
    if(max(grepl("tree_height_m", f))==0){
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
## function to clip the point cloud to a polygon and run it through the `ladderfuelsr_cbh()` function we defined above
  call_ladderfuelsr_cbh <- function(
      id
      , poly_df
      , nlas
      , my_min_vhp_n
      , my_voxel_grain_size_m
      , my_dist_btwn_bins_m
      , my_min_fuel_layer_ht_m
      , my_lad_pct_gap
      , my_lad_pct_base
      , my_num_jump_steps
      , my_min_lad_pct
      , my_frst_layer_min_ht_m
    ){
    ##################################
    # filter sf
    ##################################
    one_tree_sf <- poly_df %>% dplyr::filter(treeID==id)
    ##################################
    # clip the point cloud
    ##################################
    nlas_one_tree <- lidR::clip_roi(las = nlas, geometry = one_tree_sf) %>%
      lidR::filter_poi(!Classification %in% c(2,9,18)) %>%  ## class 2 = ground; 9 = water; 18 = noise
      lidR::add_attribute(x = id, name = "treeID")
    ##################################
    # check for points
    ##################################
    if(nrow(nlas_one_tree@data)>10){
      # safe this function
      safe_ladderfuelsr_cbh <- purrr::safely(ladderfuelsr_cbh, quiet = T)
      # # quiet that function
      quiet_ladderfuelsr_cbh <- purrr::quietly(safe_ladderfuelsr_cbh)
      # CALL ladderfuelsr_cbh
      ladderfuelsr_cbh_ans <- quiet_ladderfuelsr_cbh(
        las = nlas_one_tree
        , treeID = id
        , min_vhp_n = my_min_vhp_n
        , voxel_grain_size_m = my_voxel_grain_size_m
        , dist_btwn_bins_m = my_dist_btwn_bins_m
        , min_fuel_layer_ht_m = my_min_fuel_layer_ht_m
        , lad_pct_gap = my_lad_pct_gap
        , lad_pct_base = my_lad_pct_base
        , num_jump_steps = my_num_jump_steps
        , min_lad_pct = my_min_lad_pct
        , frst_layer_min_ht_m = my_frst_layer_min_ht_m
      )
      # get the quiet result
      ladderfuelsr_cbh_ans <- ladderfuelsr_cbh_ans$result
      # check if error
      if(is.null(ladderfuelsr_cbh_ans$error)){ # no error
        # just get the result
        ladderfuelsr_cbh_ans <- ladderfuelsr_cbh_ans$result
      }else{ # yes error
        ladderfuelsr_cbh_ans <- NULL
      }
    }else{ # < 10 points
      ladderfuelsr_cbh_ans <- NULL
    }

    # build return data
    if(is.null(ladderfuelsr_cbh_ans$cbh_metrics)){
      # blank the cbh columns
      df <- one_tree_sf %>%
        dplyr::mutate(
          cbh_maxlad_height_m = as.numeric(NA)
          , cbh_max_height_m = as.numeric(NA)
          , cbh_last_height_m = as.numeric(NA)
        )
    }else{
      df <- one_tree_sf %>%
        dplyr::mutate(
          cbh_maxlad_height_m = ladderfuelsr_cbh_ans$cbh_metrics$maxlad_Hcbh[1]
          , cbh_max_height_m = ladderfuelsr_cbh_ans$cbh_metrics$max_Hcbh[1]
          , cbh_last_height_m = ladderfuelsr_cbh_ans$cbh_metrics$last_Hcbh[1]
        )
    }
    return(df)
  }
