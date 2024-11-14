#' @title estimate CBH using `LadderFuelsR` package
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
#' @param las string -or- object. a .las|.laz file path -OR- an object of class LAS that has been height normalized.
#' @param treeID string. the LadderFuelsR process requires a treeID that uniquely identifies points within a tree, if left as NA this process will attempt to locate the `treeID` data based on an attribute in the point cloud.
#' @param min_vhp_n numeric. the minimum number of vertical height profiles (VHPs) needed to estimate a CBH.
#' @param voxel_grain_size_m numeric. horizontal resolution (suggested 1 meter for lad profiles and 10 meters for LAI maps). See `grain.size` in [leafR::lad.voxels()]
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
#' @return Returns an a list of data.frame objects that are the results of the different LadderFuelsR steps. Returns NULL if the process is unable to detect a CBH from the point cloud.
#'
#' @examples
#'  \dontrun{
#'  o <- "../data"
#'  i <- "../data/normlasdata"
#'  r <- cloud2trees::treels_stem_dbh(folder = i, outfolder = o)
#'  r %>% names()
#'  }
#' @export
#'
ladderfuelsr_cbh <- function(
  las
  , treeID = NA
  , min_vhp_n = 4
  , voxel_grain_size_m = 2
  , dist_btwn_bins_m = 1
  , min_fuel_layer_ht_m = 1
  , lad_pct_gap = 25
  , lad_pct_base = 25
  , num_jump_steps = 1
  , min_lad_pct = 10
  , frst_layer_min_ht_m = 1
) {
  # leafR does not have function reference in its coding, have to make sure this library is loaded
  library("leafR")
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
      (names(las@data) %>% stringr::str_detect("treeID") %>% max())==1
    ){
      n <- las@data$treeID %>% unique() %>% length()
      if(n>1 & is.na(treeID)){
        stop("the treeID column has more than one tree detected. set the `treeID` parameter")
      }else if(is.na(treeID)){
        # set the treeID
        treeID <- las@data$treeID %>% unique()
        # write it
        f <- las %>%
          lidR::filter_poi(treeID == treeID) %>%
          lidR::writeLAS(file = fn)
      }else{
        # write it
        f <- las %>%
          lidR::filter_poi(treeID == treeID) %>%
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
  treeID <- dplyr::coalesce(as.character(treeID), as.character(1)) # if the treeID parameter is not set, fake 1
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
      dplyr::relocate(treeID) %>%
      dplyr::filter(treeID == treeID)

  ### check if all NA or all 0, whereby no fuel gaps can be determined
  prof_na <- lad_profile %>% dplyr::filter(dplyr::coalesce(lad,0) == 0) %>% nrow()
  if( nrow(lad_profile)-prof_na <= 1 ){
    message(
      paste0(
        "no fuel gaps found. unable to quantify CBH (treeID="
        , treeID, ")."
      )
    )
    return(NULL)
  }else if(nrow(lad_profile) < min_vhp_n){
    message(
      paste0(
        nrow(lad_profile)
        , " fuel vertical height profiles found. unable to quantify CBH (treeID="
        , treeID, "). try decreasing the `min_vhp_n` parameter?"
      )
    )
    return(NULL)
  }
  else{
    ## "depurating tree lad profiles"
    ## see: https://github.com/olgaviedma/LadderFuelsR#8depurating-tree-lad-profiles
    lad_profile <- lad_profile %>%
      dplyr::mutate(lad = dplyr::coalesce(as.numeric(lad), 0.01)) %>%
      dplyr::arrange(treeID, height)
    #######################################
    ### Step 1 - `LadderFuelsR::get_gaps_fbhs`
    #######################################
      ### this function is fixed: https://github.com/olgaviedma/LadderFuelsR/pull/3
      ### LadderFuelsR::get_gaps_fbhs
      ### This function calculates gaps and fuel layers base height (FBH) as
      ### the difference in percentiles between consecutive LAD values along the vertical tree profile (VTP)

      # quiet this function
      quiet_get_gaps_fbhs <- purrr::quietly(LadderFuelsR::get_gaps_fbhs)

      gaps_fbhs <-
        # gw_get_gaps_fbhs(
        # LadderFuelsR::get_gaps_fbhs(
        quiet_get_gaps_fbhs(
          LAD_profiles = lad_profile
          , step = dist_btwn_bins_m
          , min_height = min_fuel_layer_ht_m
          , perc_gap = lad_pct_gap
          , perc_base = lad_pct_base
          , verbose = F
        )
      # just get the result
      gaps_fbhs <- gaps_fbhs$result
      # fix the columns that should be numeric
      gaps_fbhs <- gaps_fbhs %>%
        dplyr::mutate(dplyr::across(
          !tidyselect::starts_with("treeID")
          , as.numeric
        ))
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
      return(NULL)
    }else{
      ######`#################################
      ### Step 2 - `LadderFuelsR::calculate_gaps_perc`
      #######################################
        ### this function calculates the percentile value of each height
        ## LadderFuelsR::calculate_gaps_perc
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ERROR if treeID is not the first column
        # quiet this function
        quiet_calculate_gaps_perc <- purrr::quietly(LadderFuelsR::calculate_gaps_perc)
        # run it quietly
        gaps_perc <- quiet_calculate_gaps_perc(
          # LadderFuelsR::calculate_gaps_perc(
          LAD_profiles = lad_profile
          , min_height = min_fuel_layer_ht_m
        )
        # just get the result
        gaps_perc <- gaps_perc$result

      #######################################
      ### Step 3 - `LadderFuelsR::get_distance`
      #######################################
        ### calculates distances (and their heights) between fuel layers as
        ### the difference between consecutive gaps and fuel bases
        ### (the gap height always must be lower than the fuel base height).
        ## LadderFuelsR::get_distance
        metrics_distance <- LadderFuelsR::get_distance(
          gap_cbh_metrics = gaps_fbhs
          , gaps_perc = gaps_perc
          , step = dist_btwn_bins_m
          , min_height = min_fuel_layer_ht_m
          , verbose = F
        )

      #######################################
      ### Step 4 - `LadderFuelsR::get_depths`
      #######################################
        ### calculates fuels depth as the difference between gaps
        ### interleaved between fuel layers minus one step if
        ### the fuel depths are greater than one step.
        ## LadderFuelsR::get_depths
        metrics_depth <- LadderFuelsR::get_depths(
          LAD_profiles = lad_profile
          , distance_metrics = metrics_distance
          , step = dist_btwn_bins_m
          , min_height= min_fuel_layer_ht_m
          , verbose = F
        )

      #######################################
      ### Step 5 - `LadderFuelsR::get_real_fbh`
      #######################################
        ### reshapes fuel layers after removing distances equal
        ### to any number of height bin steps, keeping the first
        ### "base height" from those consecutive ones separated by such distance.

        ## LadderFuelsR::get_real_fbh
        real_fbh <- LadderFuelsR::get_real_fbh(
          depth_metrics = metrics_depth
          , step = dist_btwn_bins_m
          , number_steps = num_jump_steps
          , min_height = min_fuel_layer_ht_m
          , verbose = F
        )

      #######################################
      ### Step 6 - `LadderFuelsR::get_real_depths`
      #######################################
        ### recalculates fuel layers depth after considering
        ### distances greater than the actual height bin step.

        ## LadderFuelsR::get_real_depths
        real_depth <- LadderFuelsR::get_real_depths(
          effective_fbh = real_fbh
          , step = dist_btwn_bins_m
          , min_height = min_fuel_layer_ht_m
          , verbose = F
        )

      #######################################
      ### Step 7 - `LadderFuelsR::get_effective_gap`
      #######################################
        ### recalculates the distance between fuel layers after considering
        ### distances greater than any number of height bin steps.

        ## LadderFuelsR::get_effective_gap
        eff_gap <- LadderFuelsR::get_effective_gap(
          effective_depth = real_depth
          , number_steps = num_jump_steps
          , min_height = min_fuel_layer_ht_m
          , verbose = F
        )
      #######################################
      ### Step 8 - `LadderFuelsR::get_layers_lad`
      #######################################
        ### calculates the percentage of Leaf Area Density (LAD) within
        ### each fuel layer (first output) and removes those fuel layers
        ### with LAD percentage less than a specified threshold
        ### (default 10 the depth of the remaining ones (second output).

        ## LadderFuelsR::get_layers_lad
        layers_lad_df <- LadderFuelsR::get_layers_lad(
          LAD_profiles = lad_profile
          , effective_distances = eff_gap
          , threshold = min_lad_pct
          , step = dist_btwn_bins_m
          , min_height = min_fuel_layer_ht_m
          , verbose = F
        )
        ### idk why it is a list of 2 with the same data just the order
        ### of the `max_height` and `Hcbh1_Hdptf1` columns are switched. do you spot another difference??
        ### looking through the befuddling README it looks like the authors only keep
        ### the second data frame in the list
        if(length(layers_lad_df)>1){
          layers_lad_df <- layers_lad_df[[2]]
        }

      #######################################
      ### Step 9 - `LadderFuelsR::get_cbh_metrics`
      #######################################
        ### `LadderFuelsR::get_cbh_dist` is described in the research article but does not
        ### exist in the package or README. Looks like `LadderFuelsR::get_cbh_metrics` is there though.
        ### determines the CBH of a segmented tree using three criteria:
        ### maximum LAD percentage, maximum distance and the last distance.

        ## LadderFuelsR::get_cbh_metrics
        cbh_metrics <- LadderFuelsR::get_cbh_metrics(
          effective_LAD = layers_lad_df
          , min_height = min_fuel_layer_ht_m
          , hdepth1_height = frst_layer_min_ht_m
          , verbose = F
        )

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
    } # if all NA or all 0, whereby no fuel gaps can be determined
  } # if all NA or all 0, whereby no fuel gaps can be determined
}
