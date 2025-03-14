#' @title estimate CBH for a *single tree* using `LadderFuelsR` package
#'
#' @description
#' `ladderfuelsr_cbh()` is an all-in-one function to process height normalized .las|.laz files
#' using the functionality of the `LadderFuelsR` package.
#' The function returns a a list of data.frame objects that are the results of the different LadderFuelsR steps. Returns NULL if the process is unable to detect a CBH from the point cloud.
#' `ladderfuelsr_cbh()` outputs:
#'
#' The order of operations is:
#'
#' * Create a data frame of the 3D voxels information (xyz) with Leaf Area Density (LAD) values from las file using [leafR::lad.voxels()]
#' * Calculate the lad profile from the input lad.voxels using [leafR::lad.profile()]
#' * Calculate gaps and fuel layers base height (FBH) as the difference in percentiles between consecutive LAD values along the vertical tree profile (VTP) using [LadderFuelsR::get_gaps_fbhs()]
#' * Calculate the percentile value of each fuel layer base height using [LadderFuelsR::calculate_gaps_perc()]
#' * Calculate distances (and their heights) between fuel layers as the difference between consecutive gaps and fuel bases using [LadderFuelsR::get_distance()]
#' * Calculate fuels depth as the difference between gaps interleaved between fuel layers minus one step if the fuel depths are greater than one step using [LadderFuelsR::get_depths()]
#' * Reshape fuel layers after removing distances equal to any number of height bin steps, keeping the first "base height" from those consecutive ones separated by such distance using [LadderFuelsR::get_real_fbh()]
#' * Recalculate fuel layers depth after considering distances greater than the actual height bin step using [LadderFuelsR::get_real_depths()]
#' * Recalculate the distance between fuel layers after considering distances greater than any number of height bin steps using [LadderFuelsR::get_effective_gap()]
#' * Calculate the percentage of LAD within each fuel layer (first output) and removes those fuel layers with LAD percentage less than a specified threshold using [LadderFuelsR::get_layers_lad()]
#' * Determine the CBH of a segmented tree using three criteria: maximum LAD percentage, maximum distance and the last distance using [LadderFuelsR::get_cbh_metrics()]
#'
#' @param lad_profile_df data.frame. the return of [leafr_for_ladderfuelsr()] or a data.frame
#' that must have the columns: treeID, lad, total_pulses, height.
#' if both the `lad_profile_df` and `las` parameters are defined, the preference is the data.frame
#' @param las string -or- object. a single tree .las|.laz file path -OR- an object of class LAS that has been height normalized.
#' @param treeID numeric. the LadderFuelsR process requires a treeID that uniquely identifies points within a tree
#' , if left as NA this process will attempt to locate the `treeID`
#' data based on an attribute in the point cloud or data.frame which should be numeric
#' @param min_vhp_n numeric. the minimum number of vertical height profiles (VHPs) needed to estimate a CBH.
#' @param voxel_grain_size_m numeric. only used if `las` parameter is defined.
#' horizontal resolution (suggested 1 meter for lad profiles).
#' See `grain.size` in [leafR::lad.voxels()]
#' @param dist_btwn_bins_m numeric. value for the actual height bin step (in meters). See `step` in [LadderFuelsR::get_gaps_fbhs()]
#' @param min_fuel_layer_ht_m numeric. value for the actual minimum base height (in meters). See `min_height` in [LadderFuelsR::get_gaps_fbhs()]
#' @param lad_pct_gap numeric. value of the percentile threshold used to identify gaps (default percentile 25th). See `perc_gap` in [LadderFuelsR::get_gaps_fbhs()]
#' @param lad_pct_base numeric. value of the percentile threshold used to identify fuels layers base height (default percentile 25th). See `perc_base` in [LadderFuelsR::get_gaps_fbhs()]
#' @param num_jump_steps numeric. value for the number of height bin steps that can be jumped to reshape fuels layers. See `number_steps` in [LadderFuelsR::get_real_fbh()]
#' @param min_lad_pct numeric. value for the minimum required LAD percentage in a fuel layer. See `threshold` in [LadderFuelsR::get_layers_lad()]
#' @param frst_layer_min_ht_m numeric. value for the depth height of the first fuel layer.
#' If the first fuel layer has the maximum LAD and its depth is greater than the indicated value
#' , then this fuel layer is considered as the CBH of the tree. On the contrary, if its depth is <= the value, the CBH with maximum LAD will be the second fuel layer, although it has not the maximum LAD.
#' See `hdepth1_height` in [LadderFuelsR::get_cbh_metrics()]
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
#' @return Returns an list of data.frame objects that are the results of the different LadderFuelsR steps. Returns NULL if the process is unable to detect a CBH from the point cloud.
#'
#' @examples
#'  \dontrun{
#'  # polygon data
#'  f <- system.file(package = "cloud2trees","extdata","crowns_poly.gpkg")
#'  trees_poly <- sf::st_read(f)
#'  # simplify polygons
#'  trees_poly <- simplify_multipolygon_crowns(trees_poly)
#'  # point cloud data
#'  lf <- system.file(package = "cloud2trees","extdata","norm_las","RMNP_017_2018_normalize.las")
#'  las <- lidR::readLAS(lf)
#'  las@data %>% dplyr::glimpse()
#'  # polygon_attribute_to_las to attach treeID to las
#'  las <- polygon_attribute_to_las(las, trees_poly, force_crs = T, attribute = "treeID")
#'  las@data %>% dplyr::glimpse()
#'  # get the lad profile for each treeID
#'  lad_profile <- leafr_for_ladderfuelsr(
#'      las
#'      , voxel_grain_size_m = 1
#'      , k = 1
#'      , group_treeID = T
#'      , relative = F
#'    )
#'  dplyr::glimpse(lad_profile)
#'  # extract the CBH using ladderfuelsr_cbh()
#'  # before we extract the CBH using ladderfuelsr_cbh(), treeID has to be numeric
#'  lad_profile <- lad_profile %>%
#'    dplyr::mutate(
#'      treeID_backup = treeID, treeID = as.factor(treeID)
#'    )
#'  # for one tree
#'  ladderfuelsr_cbh(
#'      lad_profile_df = lad_profile
#'      , treeID = lad_profile$treeID[1]
#'    ) %>%
#'    purrr::pluck("cbh_metrics") %>%
#'    dplyr::glimpse()
#'  # we can map over multiple trees
#'  cbhs <- lad_profile$treeID %>%
#'    unique() %>%
#'    .[1:22] %>%
#'    purrr::map(\(x) ladderfuelsr_cbh(
#'        lad_profile_df = lad_profile
#'        , treeID = x
#'      ) %>%
#'      purrr::pluck("cbh_metrics")
#'    ) %>%
#'    dplyr::bind_rows()
#'  dplyr::glimpse(cbhs)
#'  ggplot2::ggplot(data = cbhs, mapping = ggplot2::aes(x=last_Hcbh)) +
#'    ggplot2::geom_density()
#'  }
#' @export
#'
#'
ladderfuelsr_cbh <- function(
  lad_profile_df = NULL
  , las = NULL
  , treeID = NA
  , min_vhp_n = 4
  , voxel_grain_size_m = 1
  , dist_btwn_bins_m = 1
  , min_fuel_layer_ht_m = 1
  , lad_pct_gap = 25
  , lad_pct_base = 25
  , num_jump_steps = 1
  , min_lad_pct = 10
  , frst_layer_min_ht_m = 1
) {
  # LadderFuelsR required
    if(!requireNamespace("LadderFuelsR", quietly = TRUE)) {
      stop(paste0(
        "Package \"LadderFuelsR\" must be installed to use this function."
        , "\n"
        , "try `pak::pak(\"olgaviedma/LadderFuelsR\", upgrade = TRUE)`"
      ))
    }
  # set up empty return
    empty_return <- list(
      gaps_fbhs = NULL
      , lad_profile = NULL
      , gaps_perc = NULL
      , metrics_distance = NULL
      , metrics_depth = NULL
      , real_fbh = NULL
      , real_depth = NULL
      , eff_gap = NULL
      , layers_lad_df = NULL
      , cbh_metrics = NULL
    )
  ###########################################################################
  ###########################################################################
  # setup:
  #   A) if lad_profile_df is a data.frame:
  #     - filter for the treeID if defined,
  #       otherwise if more than treeID detected will throw error
  #     - checks required columns
  #     - "depurating tree lad profiles"
  #     - the lad_profile data.frame goes to the LadderFuelsR steps
  #   b) if las is a character or LAS
  #     - checks if can read .las/.laz from file or is LAS, stops if not
  #     - filter for the treeID if defined,
  #       otherwise if more than treeID detected will throw error
  #     - uses leafR package to get lad profiles from the point cloud
  #     - "depurating tree lad profiles"
  #     - the lad_profile data.frame goes to the LadderFuelsR steps
  ###########################################################################
  ###########################################################################
  if(inherits(lad_profile_df, "data.frame")){
    # check treeID
    if(
      !(names(lad_profile_df) %>% stringr::str_equal("treeID") %>% any())
    ){
      # set the treeID
      lad_profile_df <- lad_profile_df %>%
        dplyr::mutate(
          # if the treeID parameter is not set, fake 1
          treeID = dplyr::coalesce(
              as.numeric(treeID)
              , round(as.numeric(Sys.time()) + round(runif(1)*100000))
            ) %>% 
            as_character_safe() %>% 
            as.factor()
        )
    }else if( # filter for the treeID
      (names(lad_profile_df) %>% stringr::str_equal("treeID") %>% any())
      && !is.na(treeID) && !is.null(treeID)
    ){
      id_temp <- as.numeric(treeID)[1]
      lad_profile_df <- lad_profile_df %>%
        dplyr::filter(as.numeric(treeID)==id_temp)
      if(nrow(lad_profile_df)<1){
        stop("the treeID was not detected in the data. double-check the `treeID` parameter or leave as NA")
      }
    }

    #check the required columns which should have been
    # exported from leafr_for_ladderfuelsr()
    req_cols <- c(
      "treeID"
      , "lad"
      , "total_pulses"
      , "height"
    )
    # check_df_cols_all_missing() in utils_biomass.r...will stop program if any missing
    check_df_cols_all_missing(
      lad_profile_df
      , col_names = req_cols %>% unique()
      , check_vals_missing = F
    )
    # filter for trees where CBH process implemented via LadderFuelsR will be successful
    ## "depurating tree lad profiles"
    ## see: https://github.com/olgaviedma/LadderFuelsR#8depurating-tree-lad-profiles
    lad_profile <-
      lad_profile_df %>%
      dplyr::group_by(treeID) %>%
      # count the number of vertical height profiles
      dplyr::mutate(
        vhp = sum( ifelse(dplyr::coalesce(as.numeric(lad),0)>0, 1, 0) )
      ) %>%
      dplyr::filter(
        # the LadderFuelsR documentation filters for .las files with more than 10 points
        as.numeric(total_pulses) >= 6
        # no fuel gaps can be determined if < 1 vhps with >0 lad
        & vhp > 0
        & dplyr::n() >= 2 ## have to have at least 2 vhps in addition to at least one vhp>0
        & dplyr::n() >= dplyr::coalesce(as.numeric(min_vhp_n),2) #user defined minimum
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        lad = dplyr::coalesce(as.numeric(lad), 0.01)
        , height = as.numeric(height)
      ) %>%
      ## !!!!! not only does the treeID column have to exist...it has to be the first column
      ## !!!!! what if treeID is not numeric? idk???
      dplyr::select(treeID, height, lad) %>%
      dplyr::arrange(treeID, height)

    ## check rows and return null if none
    if(nrow(lad_profile)<1){return(empty_return)}

    ## check for multiple treeIDs
    if( length(unique(lad_profile$treeID))>1 ){
      stop("lad profile data contains mulitple trees, set the `treeID` parameter or filter data for a single tree first")
    }

  }else if(
    inherits(las, "character") ||
    inherits(las, "LAS")
  ){ # checks if passed a las file instead of data.frame
    # leafR required
      if(!requireNamespace("leafR", quietly = TRUE)) {
        stop(paste0(
          "Package \"leafR\" must be installed to use this function."
          , "\n"
          , "try `pak::pak(\"DRAAlmeida/leafR\", upgrade = TRUE)`"
        ))
      }

      # check if string to las/laz file
      if(inherits(las, "character")){
        if(!stringr::str_ends(las, ".*\\.(laz|las)$")){
          stop("must pass a .las|.laz file path -OR- an object of class LAS to the `las` parameter")
        }
        # set the file path
        f <- normalizePath(las)
      }else if(inherits(las, "LAS")){
        # have to write the las to a tempfile
        fn <- paste0(tempdir(), "/temp.las")
        # check if has a treeID
        if(
          (names(las@data) %>% stringr::str_detect("treeID") %>% any())
        ){
          # filter for the treeID
          if(!is.na(treeID) && !is.null(treeID)){
            id_temp <- as.numeric(treeID)[1]
            las <- las %>% lidR::filter_poi(treeID == id_temp)
            if(lidR::is.empty(las)){
              stop("the treeID was not detected in the data. double-check the `treeID` parameter or leave as NA")
            }
          }
          # check for multiple treeID
          n <- las@data$treeID %>% unique() %>% length()
          if(n>1 & is.na(treeID)){
            stop("the treeID column has more than one tree detected. set the `treeID` parameter")
          }else if(is.na(treeID)){
            # set the treeID
            id_temp <- las@data$treeID %>% unique()
            treeID <- id_temp
            # write it
            f <- las %>%
              lidR::filter_poi(treeID == id_temp) %>%
              lidR::writeLAS(file = fn)
          }else{
            id_temp <- as.numeric(treeID)[1]
            # write it
            f <- las %>%
              lidR::filter_poi(treeID == id_temp) %>%
              lidR::writeLAS(file = fn)
          }
        }else{
          # write it
          f <- las %>% lidR::writeLAS(file = fn)
        }
      }else{
        stop("must pass a .las|.laz file path -OR- an object of class LAS to the `las` parameter")
      }

      # check the treeID
      treeID <- dplyr::coalesce(
        as.numeric(treeID)[1]
        , round(as.numeric(Sys.time()) + round(runif(1)*100000))
      ) %>% 
      as_character_safe() %>% 
      as.factor() # if the treeID parameter is not set, fake 1
      #######################################
      ### Step 0 - `leafR` steps
      #######################################
        # 1) `leafR::lad.voxels()` - use normalized las file to create
            # a data frame of the 3D voxels information (xyz) with Leaf Area Density values
        # 2) `leafR::lad.profile()` - calculate the lad profile from
            # the input lad.voxels (step 1)
        # 3) ensure that the data frame returned from `leafR::lad.profile()`
            # has a column named `treeID` which uniquely identifies individual trees.
            # also, that column has to be the first column (bad practice by the authors)

        ## leafR::lad.voxels
        lad_voxels <- leafR::lad.voxels(normlas.file = f, grain.size = voxel_grain_size_m)
        ## leafR::lad.profile
        lad_profile <- leafR::lad.profile(lad_voxels, relative = F)
        ## add treeID column that is required by the package, though it's never stated
        lad_profile <- lad_profile %>%
          dplyr::mutate(
            treeID = treeID %>% factor()
          ) %>%
          ## !!!!! not only does the treeID column have to exist...it has to be the first column
          dplyr::relocate(treeID)

      ### check if all NA or all 0, whereby no fuel gaps can be determined
      prof_na <- lad_profile %>% dplyr::filter(dplyr::coalesce(lad,0) == 0) %>% nrow()
      if( nrow(lad_profile)-prof_na <= 1 ){
        message(
          paste0(
            "no fuel gaps found. unable to quantify CBH (treeID="
            , treeID, ")."
          )
        )
        empty_return$lad_profile <- lad_profile
        return(empty_return)
      }else if(nrow(lad_profile) < min_vhp_n){
        message(
          paste0(
            nrow(lad_profile)
            , " fuel vertical height profiles found. unable to quantify CBH (treeID="
            , treeID, "). try decreasing the `min_vhp_n` parameter?"
          )
        )
        empty_return$lad_profile <- lad_profile
        return(empty_return)
      }
      else{
        ## "depurating tree lad profiles"
        ## see: https://github.com/olgaviedma/LadderFuelsR#8depurating-tree-lad-profiles
        lad_profile <- lad_profile %>%
          dplyr::mutate(lad = dplyr::coalesce(as.numeric(lad), 0.01)) %>%
          dplyr::arrange(treeID, height)
      }

  }else{ # check if passed las or data.frame
    stop("must define `lad_profile_df` or `las` parameter")
  }

  ###########################################################################
  ###########################################################################
  # the lad_profile data.frame goes to the LadderFuelsR steps
  ###########################################################################
  ###########################################################################

    #######################################
    ### Step 1 - `LadderFuelsR::get_gaps_fbhs`
    #######################################
      ### this function is fixed: https://github.com/olgaviedma/LadderFuelsR/pull/3
      ### LadderFuelsR::get_gaps_fbhs
      ### This function calculates gaps and fuel layers base height (FBH) as
      ### the difference in percentiles between consecutive LAD values along the vertical tree profile (VTP)
      ## LadderFuelsR::get_gaps_fbhs
      # safe this function
      safe_get_gaps_fbhs <- purrr::safely(LadderFuelsR::get_gaps_fbhs)

      gaps_fbhs <-
        safe_get_gaps_fbhs(
          LAD_profiles = lad_profile
          , step = dist_btwn_bins_m
          , min_height = min_fuel_layer_ht_m
          , perc_gap = lad_pct_gap
          , perc_base = lad_pct_base
          , verbose = F
        )
      # just get the result
      gaps_fbhs <- gaps_fbhs$result

      if(dplyr::coalesce(nrow(gaps_fbhs),0)<1){
        empty_return$lad_profile <- lad_profile
        return(empty_return)
      }
      # fix the columns that should be numeric
        gaps_fbhs <- gaps_fbhs %>%
          dplyr::mutate(dplyr::across(
            !tidyselect::starts_with("treeID")
            , as.numeric
          )) %>%
          dplyr::relocate(treeID)
      ### check for all NA or 0
      gaps_na <- gaps_fbhs %>%
        dplyr::filter(
          dplyr::if_all(
            .cols = -tidyselect::starts_with("treeID")
            , .fns = ~ dplyr::coalesce(.x, 0) == 0
          )
        ) %>%
        nrow()
      if(gaps_na>0){
        message(
          paste0(
            "no fuel gaps found. unable to quantify CBH (treeID="
            , treeID, ")."
          )
        )

        empty_return$lad_profile <- lad_profile
        return(empty_return)
      }

    #######################################
    ### Step 2 - `LadderFuelsR::calculate_gaps_perc`
    #######################################
      ### this function calculates the percentile value of each height
      ## LadderFuelsR::calculate_gaps_perc
      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ERROR if treeID is not the first column
      # safe this function
      safe_calculate_gaps_perc <- purrr::safely(LadderFuelsR::calculate_gaps_perc)
      # run it safely
      gaps_perc <- safe_calculate_gaps_perc(
        # LadderFuelsR::calculate_gaps_perc(
        LAD_profiles = lad_profile
        , min_height = min_fuel_layer_ht_m
      )
      # just get the result
      gaps_perc <- gaps_perc$result
      # stop if nothing found
      if(dplyr::coalesce(nrow(gaps_perc),0)<1){
        empty_return$lad_profile <- lad_profile
        empty_return$gaps_fbhs <- gaps_fbhs
        return(empty_return)
      }
      # fix the columns that should be numeric
      gaps_perc <- gaps_perc %>%
        dplyr::mutate(dplyr::across(
          !tidyselect::starts_with("treeID")
          , as.numeric
        )) %>%
        dplyr::relocate(treeID)
    #######################################
    ### Step 3 - `LadderFuelsR::get_distance`
    #######################################
      ### calculates distances (and their heights) between fuel layers as
      ### the difference between consecutive gaps and fuel bases
      ### (the gap height always must be lower than the fuel base height).
      ## LadderFuelsR::get_distance

      # safe this function
      safe_get_distance <- purrr::safely(LadderFuelsR::get_distance)
      # run it safely
      metrics_distance <- safe_get_distance(
        gap_cbh_metrics = gaps_fbhs
        , gaps_perc = gaps_perc
        , step = dist_btwn_bins_m
        , min_height = min_fuel_layer_ht_m
        , verbose = F
      )
      # just get the result
      metrics_distance <- metrics_distance$result
      # stop if nothing found
      if(dplyr::coalesce(nrow(metrics_distance),0)<1){
        empty_return$lad_profile <- lad_profile
        empty_return$gaps_fbhs <- gaps_fbhs
        empty_return$gaps_perc <- gaps_perc
        return(empty_return)
      }
      # fix the columns that should be numeric
      metrics_distance <- metrics_distance %>%
        dplyr::mutate(dplyr::across(
          !tidyselect::starts_with("treeID")
          , as.numeric
        )) %>%
        dplyr::relocate(treeID)
    #######################################
    ### Step 4 - `LadderFuelsR::get_depths`
    #######################################
      ### calculates fuels depth as the difference between gaps
      ### interleaved between fuel layers minus one step if
      ### the fuel depths are greater than one step.
      ## LadderFuelsR::get_depths

      # safe this function
      ### there is an error that needs to be fixed
      ### directly in LadderFuelsR::get_depths:
        ### Error in data.frame(..., check.names = FALSE) :
        ### arguments imply differing number of rows: 1, 0
      safe_get_depths <- purrr::safely(LadderFuelsR::get_depths)

      metrics_depth <- safe_get_depths(
        LAD_profiles = lad_profile
        , distance_metrics = metrics_distance
        , step = dist_btwn_bins_m
        , min_height= min_fuel_layer_ht_m
        , verbose = F
      )
      # just get the result
      metrics_depth <- metrics_depth$result
      # stop if nothing found
      if(dplyr::coalesce(nrow(metrics_depth),0)<1){
        empty_return$lad_profile <- lad_profile
        empty_return$gaps_fbhs <- gaps_fbhs
        empty_return$gaps_perc <- gaps_perc
        empty_return$metrics_distance <- metrics_distance
        return(empty_return)
      }
      # fix the columns that should be numeric
      metrics_depth <- metrics_depth %>%
        dplyr::mutate(dplyr::across(
          !tidyselect::starts_with("treeID")
          , as.numeric
        )) %>%
        dplyr::relocate(treeID)
    #######################################
    ### Step 5 - `LadderFuelsR::get_real_fbh`
    #######################################
      ### reshapes fuel layers after removing distances equal
      ### to any number of height bin steps, keeping the first
      ### "base height" from those consecutive ones separated by such distance.
      ## LadderFuelsR::get_real_fbh

      # safe this function
      safe_get_real_fbh <- purrr::safely(LadderFuelsR::get_real_fbh)
      # run it safely
      real_fbh <- safe_get_real_fbh(
        depth_metrics = metrics_depth
        , step = dist_btwn_bins_m
        , number_steps = num_jump_steps
        , min_height = min_fuel_layer_ht_m
        , verbose = F
      )
      # just get the result
      real_fbh <- real_fbh$result
      # stop if nothing found
      if(dplyr::coalesce(nrow(real_fbh),0)<1){
        empty_return$lad_profile <- lad_profile
        empty_return$gaps_fbhs <- gaps_fbhs
        empty_return$gaps_perc <- gaps_perc
        empty_return$metrics_distance <- metrics_distance
        empty_return$metrics_depth <- metrics_depth
        return(empty_return)
      }
      # fix the columns that should be numeric
      real_fbh <- real_fbh %>%
        dplyr::mutate(dplyr::across(
          !tidyselect::starts_with("treeID")
          , as.numeric
        )) %>%
        dplyr::relocate(treeID)

    #######################################
    ### Step 6 - `LadderFuelsR::get_real_depths`
    #######################################
      ### recalculates fuel layers depth after considering
      ### distances greater than the actual height bin step.
      ## LadderFuelsR::get_real_depths

      # safe this function
      safe_get_real_depths <- purrr::safely(LadderFuelsR::get_real_depths)
      # run it safely
      real_depth <- safe_get_real_depths(
        effective_fbh = real_fbh
        , step = dist_btwn_bins_m
        , min_height = min_fuel_layer_ht_m
        , verbose = F
      )
      # just get the result
      real_depth <- real_depth$result
      # stop if nothing found
      if(dplyr::coalesce(nrow(real_depth),0)<1){
        empty_return$lad_profile <- lad_profile
        empty_return$gaps_fbhs <- gaps_fbhs
        empty_return$gaps_perc <- gaps_perc
        empty_return$metrics_distance <- metrics_distance
        empty_return$metrics_depth <- metrics_depth
        empty_return$real_fbh <- real_fbh
        return(empty_return)
      }
      # fix the columns that should be numeric
      real_depth <- real_depth %>%
        dplyr::mutate(dplyr::across(
          !tidyselect::starts_with("treeID")
          , as.numeric
        )) %>%
        dplyr::relocate(treeID)

    #######################################
    ### Step 7 - `LadderFuelsR::get_effective_gap`
    #######################################
      ### recalculates the distance between fuel layers after considering
      ### distances greater than any number of height bin steps.
      ## LadderFuelsR::get_effective_gap

      # safe this function
      safe_get_effective_gap <- purrr::safely(LadderFuelsR::get_effective_gap)
      # run it safely
      eff_gap <- safe_get_effective_gap(
        effective_depth = real_depth
        , number_steps = num_jump_steps
        , min_height = min_fuel_layer_ht_m
        , verbose = F
      )
      # just get the result
      eff_gap <- eff_gap$result
      # stop if nothing found
      if(dplyr::coalesce(nrow(eff_gap),0)<1){
        empty_return$lad_profile <- lad_profile
        empty_return$gaps_fbhs <- gaps_fbhs
        empty_return$gaps_perc <- gaps_perc
        empty_return$metrics_distance <- metrics_distance
        empty_return$metrics_depth <- metrics_depth
        empty_return$real_fbh <- real_fbh
        empty_return$real_depth <- real_depth
        return(empty_return)
      }
      # fix the columns that should be numeric
      eff_gap <- eff_gap %>%
        dplyr::mutate(dplyr::across(
          !tidyselect::starts_with("treeID")
          , as.numeric
        )) %>%
        dplyr::relocate(treeID)

    #######################################
    ### Step 8 - `LadderFuelsR::get_layers_lad`
    #######################################
      ### calculates the percentage of Leaf Area Density (LAD) within
      ### each fuel layer (first output) and removes those fuel layers
      ### with LAD percentage less than a specified threshold
      ### (default 10 the depth of the remaining ones (second output).
      ## LadderFuelsR::get_layers_lad

      # safe this function
      safe_get_layers_lad <- purrr::safely(LadderFuelsR::get_layers_lad)
      # run it safely
      layers_lad_df <- safe_get_layers_lad(
        LAD_profiles = lad_profile
        , effective_distances = eff_gap
        , threshold = min_lad_pct
        , step = dist_btwn_bins_m
        , min_height = min_fuel_layer_ht_m
        , verbose = F
      )
      # just get the result
      layers_lad_df <- layers_lad_df$result
      ### idk why it is a list of 2 with the same data just the order
      ### of the `max_height` and `Hcbh1_Hdptf1` columns are switched. do you spot another difference??
      ### looking through the befuddling README it looks like the authors only keep
      ### the second data frame in the list
      if(length(layers_lad_df)>1){
        layers_lad_df <- layers_lad_df[[2]]
      }
      # stop if nothing found
      if(dplyr::coalesce(nrow(layers_lad_df),0)<1){
        empty_return$lad_profile <- lad_profile
        empty_return$gaps_fbhs <- gaps_fbhs
        empty_return$gaps_perc <- gaps_perc
        empty_return$metrics_distance <- metrics_distance
        empty_return$metrics_depth <- metrics_depth
        empty_return$real_fbh <- real_fbh
        empty_return$real_depth <- real_depth
        empty_return$eff_gap <- eff_gap
        return(empty_return)
      }
      # fix the columns that should be numeric
      layers_lad_df <- layers_lad_df %>%
        dplyr::mutate(dplyr::across(
          !tidyselect::starts_with("treeID")
          , as.numeric
        )) %>%
        dplyr::relocate(treeID)

    #######################################
    ### Step 9 - `LadderFuelsR::get_cbh_metrics`
    #######################################
      ### `LadderFuelsR::get_cbh_dist` is described in the research article but does not
      ### exist in the package or README. Looks like `LadderFuelsR::get_cbh_metrics` is there though.
      ### determines the CBH of a segmented tree using three criteria:
      ### maximum LAD percentage, maximum distance and the last distance.
      ## LadderFuelsR::get_cbh_metrics

      # safe this function
      safe_get_cbh_metrics <- purrr::safely(LadderFuelsR::get_cbh_metrics)
      # run it safely
      cbh_metrics <- safe_get_cbh_metrics(
        effective_LAD = layers_lad_df
        , min_height = min_fuel_layer_ht_m
        , hdepth1_height = frst_layer_min_ht_m
        , verbose = F
      )
      # just get the result
      cbh_metrics <- cbh_metrics$result
      # stop if nothing found
      if(dplyr::coalesce(nrow(cbh_metrics),0)<1){
        empty_return$lad_profile <- lad_profile
        empty_return$gaps_fbhs <- gaps_fbhs
        empty_return$gaps_perc <- gaps_perc
        empty_return$metrics_distance <- metrics_distance
        empty_return$metrics_depth <- metrics_depth
        empty_return$real_fbh <- real_fbh
        empty_return$real_depth <- real_depth
        empty_return$layers_lad_df <- layers_lad_df
        return(empty_return)
      }
      # fix the columns that should be numeric
      cbh_metrics <- cbh_metrics %>%
        dplyr::mutate(dplyr::across(
          !tidyselect::starts_with("treeID")
          , as.numeric
        )) %>%
        dplyr::relocate(treeID)

    # return
      return(list(
        gaps_fbhs = gaps_fbhs
        , lad_profile = lad_profile
        , gaps_perc = gaps_perc
        , metrics_distance = metrics_distance
        , metrics_depth = metrics_depth
        , real_fbh = real_fbh
        , real_depth = real_depth
        , eff_gap = eff_gap
        , layers_lad_df = layers_lad_df
        , cbh_metrics = cbh_metrics
      ))

}
#############################################################
# leafR does not have function reference in its coding and calls this nefarious "pointsByZSlice"
pointsByZSlice <<- leafR::pointsByZSlice
# # see: https://github.com/DRAAlmeida/leafR/blob/fd1456b9692ba7b5d5ba94ae88216046c8ec186f/R/main.R#L11
## sets this function as global so that it works in leafR commands
# pointsByZSlice <<- function(Z, maxZ){
#     heightSlices = as.integer(Z) # Round down
#     zSlice = data.table::data.table(Z=Z, heightSlices=heightSlices) # Create a data.table (Z, slices))
#     sliceCount = stats::aggregate(list(V1=Z), list(heightSlices=heightSlices), length) # Count number of returns by slice
#
#     ##############################################
#     # Add columns to equalize number of columns
#     ##############################################
#     colRange = 0:maxZ
#     addToList = setdiff(colRange, sliceCount$heightSlices)
#     n = length(addToList)
#     if (n > 0) {
#       bindDt = data.frame(heightSlices = addToList, V1=integer(n))
#       sliceCount = rbind(sliceCount, bindDt)
#       # Order by height
#       sliceCount = sliceCount[order(sliceCount$heightSlices),]
#     }
#
#     colNames = as.character(sliceCount$heightSlices)
#     colNames[1] = "ground_0_1m"
#     colNames[-1] = paste0("pulses_", colNames[-1], "_", sliceCount$heightSlices[-1]+1, "m")
#     metrics = list()
#     metrics[colNames] = sliceCount$V1
#
#     return(metrics)
#
#   } #end function pointsByZSlice
