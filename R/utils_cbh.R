#' @title Estimate CBH using tree crown polygons and normalized point cloud data
#'
#' @description
#' `trees_cbh_sf()` uses the input tree crown polygons (e.g. as exported by [raster2trees()]) with the columns
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
#' Or the path to a single spatial polygon file.
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
#' @param trees_sample data.frame. provide your own tree sample list such as one generated from `sample_trees_flist()` that includes the `treeID` column.
#' If provided, the tree_sample_n,tree_sample_prop will be ignored
#' @param ofile character or logical. if a character value is provided the output will be written to the disk as a csv at the location provided.
#' If set to TRUE and a file path was used as the input for `trees_poly`, then a csv file will be written to the same location with the same name prefixed with "cbh_".
#' Leave as NA to return a data.frame of the trees from tree list from `trees_poly` with CBH values added
#'
#' @keywords internal
#'
trees_cbh_sf <- function(
  trees_poly
  , norm_las = NULL
  , tree_sample_n = NA
  , tree_sample_prop = NA
  , which_cbh = "lowest"
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
  , trees_sample = NA # answer from sample_trees_flist()
  , ofile = NA
){
  force_cbh_lte_ht <- T
  ##################################
  # check sample proportion
  ##################################
    check_sample_vals_ans <- check_sample_vals(
      tree_sample_n = tree_sample_n
      , tree_sample_prop = tree_sample_prop
      , def_tree_sample_n = 333
    )
    tree_sample_n <- check_sample_vals_ans$tree_sample_n
    tree_sample_prop <- check_sample_vals_ans$tree_sample_prop
  ##################################
  # check which cbh
  ##################################
    which_cbh <- check_which_cbh(which_cbh = which_cbh)
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
  # check if passed a file name or sf
  ##################################
    # message
    sf_msg <- paste0(
      "`trees_poly` data must be an object of class `sf` with only POLYGON type."
      , "\nProvide an `sf` object and see `sf::st_geometry_type()`."
    )
    # blank parameter to save backup of trees_poly if it's a character
    fn_for_ofile <- NULL
    if(
      !inherits(trees_poly, "sf")
      && inherits(trees_poly, "character")
      && file.exists(trees_poly[1]) # its a findable file
      && !dir.exists(trees_poly[1]) # its not a directory
    ){
      # save the file name
      fn_for_ofile <- normalizePath(trees_poly[1])
      # read crown
      trees_poly <- sf::st_read(dsn = fn_for_ofile, quiet = T)
      # check
      if(nrow(trees_poly)==0){return(sf_msg)}
    }else if(inherits(trees_poly, "character")){
      stop(paste0(
        "the character value provided in `trees_poly` is not a spatial file at \n    "
        , trees_poly
        , "\n   use trees_cbh_flist() instead?"
      ))
    }
  ##################################
  # ensure that treeID data exists
  ##################################
  f <- trees_poly %>% names() %>% dplyr::coalesce("")
  # check_trees_poly will throw error if fails any checks
  check_trees_poly_ans <- check_trees_poly(trees_poly, orig_fnm = fn_for_ofile)
  id_class <- class(trees_poly$treeID)[1]

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
      inherits(trees_sample, "data.frame")
      && (names(trees_sample) %>% stringr::str_equal("treeID") %>% any())
    ){
      samp_trees <- trees_poly %>%
        # this keeps our sample only
        dplyr::inner_join(trees_sample, by = "treeID")
      if(nrow(samp_trees)==0){return(NULL)}
    }else if(
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

  # check if we need to split for massive tree crown data
  if(nrow(simp_trees_poly)>500e3){
    # for data with so many crowns I've encountered the error:
    ### !!! Error in getGlobalsAndPackages(expr, envir = envir, tweak = tweakExpression, :
      ### !!! The total size of the xx globals exported for future expression ...
      ### !!! is xxx GiB.. This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize').
    # break up data
    simp_trees_poly <-
      simp_trees_poly %>%
      # arrange by x,y
      dplyr::left_join(
        simp_trees_poly %>%
          sf::st_centroid() %>%
          dplyr::mutate(
            x_xxx = sf::st_coordinates(.)[,1]
            , y_xxx = sf::st_coordinates(.)[,2]
          ) %>%
          sf::st_drop_geometry() %>%
          dplyr::select(treeID, x_xxx, y_xxx)
        , by = "treeID"
      ) %>%
      dplyr::arrange(x_xxx,y_xxx) %>%
      # groups of 250k....or larger
      dplyr::mutate(grp = ceiling(dplyr::row_number()/500e3)) %>%
      dplyr::select(-c(x_xxx,y_xxx))

    # apply it
    output_temp <- simp_trees_poly$grp %>%
      unique() %>%
      purrr::map(function(x){
        # rename output
        lidR::opt_output_files(nlas_ctg) <- paste0(tempdir(), "/{*}_treed_",x)
        # run it
        lidR::catalog_apply(
          ctg = nlas_ctg
          , FUN = ctg_leafr_for_ladderfuelsr
          , .options = list(automerge = TRUE)
          # ctg_calc_tree_cbh options
          , poly_df = simp_trees_poly %>% dplyr::filter(grp==x)
          , force_crs = force_same_crs
          , voxel_grain_size_m = voxel_grain_size_m
        )
      }) %>%
      unlist()
  }else{
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
  }
  ##################################
  # read and clean output files for ladderfuelsr_cbh
  ##################################
  lad_profile <- output_ctg_for_ladderfuelsr_cbh(
    output = output_temp
    , id_class = id_class
  )
  ##################################
  # do ladderfuelsr_cbh
  ##################################
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
  ##################################
  # check the cbh data we got and use the cbh selected
  ##################################
  cbh_df <- clean_cbh_df(
    cbh_df = cbh_df
    , trees_poly = trees_poly
    , lad_profile = lad_profile
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
      if(
        new_voxel_grain_size_m!=voxel_grain_size_m
      ){
        # check if we need to split for massive tree crown data
        if(nrow(simp_trees_poly)>500e3){
          # for data with so many crowns I've encountered the error:
          ### !!! Error in getGlobalsAndPackages(expr, envir = envir, tweak = tweakExpression, :
            ### !!! The total size of the xx globals exported for future expression ...
            ### !!! is xxx GiB.. This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize').
          # apply it
          output_temp <- simp_trees_poly$grp %>%
            unique() %>%
            purrr::map(function(x){
              # rename output
              lidR::opt_output_files(nlas_ctg) <- paste0(tempdir(), "/{*}_treed_",x)
              # run it
              output_temp <- lidR::catalog_apply(
                ctg = nlas_ctg
                , FUN = ctg_leafr_for_ladderfuelsr
                , .options = list(automerge = TRUE)
                # ctg_calc_tree_cbh options
                , poly_df = simp_trees_poly %>% dplyr::filter(grp==x)
                , force_crs = force_same_crs
                , voxel_grain_size_m = new_voxel_grain_size_m
              )
            }) %>%
            unlist()
        }else{
          # apply it
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
      }

      # did we change anything?
      if(
        min_fuel_layer_ht_m != new_min_fuel_layer_ht_m
        || dist_btwn_bins_m != new_dist_btwn_bins_m
        || min_vhp_n != new_min_vhp_n
        || voxel_grain_size_m != new_voxel_grain_size_m
        || min_lad_pct != new_min_lad_pct
      ){
        ##################################
        # read and clean output files for ladderfuelsr_cbh
        ##################################
        lad_profile <- output_ctg_for_ladderfuelsr_cbh(
          output = output_temp
          , id_class = id_class
        )
        ##################################
        # do ladderfuelsr_cbh
        ##################################
        # dplyr::glimpse(lad_profile)
        # extract the CBH using ladderfuelsr_cbh()
        # we can map over multiple trees

        cbh_df <- lad_profile$treeID %>%
          unique() %>%
          purrr::map(\(x)
              quiet_ladderfuelsr_cbh(
              # ladderfuelsr_cbh(
                lad_profile_df = lad_profile
                , treeID = x
                , min_vhp_n = new_min_vhp_n
                , voxel_grain_size_m = new_voxel_grain_size_m
                , dist_btwn_bins_m = new_dist_btwn_bins_m
                , min_fuel_layer_ht_m = new_min_fuel_layer_ht_m
                , lad_pct_gap = lad_pct_gap
                , lad_pct_base = lad_pct_base
                , num_jump_steps = num_jump_steps
                , min_lad_pct = new_min_lad_pct
                , frst_layer_min_ht_m = frst_layer_min_ht_m
              ) %>%
              purrr::pluck("result") %>% ## b/c purrr::quietly
              purrr::pluck("cbh_metrics")
            , .progress = "extracting CBH"
          ) %>%
          dplyr::bind_rows()
        # dplyr::glimpse(cbh_df)
        ##################################
        # check the cbh data we got and use the cbh selected
        ##################################
        cbh_df <- clean_cbh_df(
          cbh_df = cbh_df
          , trees_poly = trees_poly
          , lad_profile = lad_profile
          , force_cbh_lte_ht = force_cbh_lte_ht
          , which_cbh = which_cbh
        )

        # ensure that there are enough data to estimate
        n_cbh <- nrow(cbh_df)
      }
    }

  ##################################
  # check if write it and return
  ##################################
  if(
    !is.na(ofile)
    && !is.null(ofile)
    && ofile==T
    && !is.null(fn_for_ofile)
    && !is.na(fn_for_ofile)
  ){
    fnm <- paste0(
      "cbh_"
      , fn_for_ofile %>%
        basename() %>%
        # remove the file extension
        stringr::str_replace(pattern = "(.*)\\..*$", replacement = "\\1")
      , ".csv"
    )
    fpth <- file.path(
      dirname(fn_for_ofile)
      , fnm
    )
    # write it
    cbh_df %>%
      sf::st_drop_geometry() %>%
      write.csv(file = fpth, row.names = F, append = F)
    # return
    return(fpth)
  }else if(
    !is.na(ofile)
    && inherits(ofile,"character")
  ){
    if(
      # you sent me a file name
      stringr::str_ends(ofile, ".*\\.csv$")
    ){
      # write it
      cbh_df %>%
        sf::st_drop_geometry() %>%
        write.csv(file = ofile, row.names = F, append = F)
      # return
      return(ofile)
    }else if(
      # you sent me a directory
      dir.exists( normalizePath(ofile) )
    ){
      fpth <- file.path(
        normalizePath(ofile)
        , "cbh_tree_list.csv"
      )
      # write it
      cbh_df %>%
        sf::st_drop_geometry() %>%
        write.csv(file = fpth, row.names = F, append = F)
      # return
      return(fpth)
    }else{ # you sent me a filename without .csv
      fpth <- paste0(
        normalizePath(ofile)
        , ".csv"
      )
      # write it
      cbh_df %>%
        sf::st_drop_geometry() %>%
        write.csv(file = fpth, row.names = F, append = F)
      # return
      return(fpth)
    }
  }else{
    # return
    return(cbh_df)
  }
}

################################################################################################################################################
# function to apply trees_cbh_sf to a directory where cloud2trees::cloud2trees() output was written
# or to a list of crown files
# this function pipeline was developed to avoid memory issues with xxl tree lists
################################################################################################################################################
trees_cbh_flist <- function(
  # flist == a list of crown files or the directory where final_detected_tree_tops* and final_detected_crowns* files
  # from cloud2trees::cloud2trees() or cloud2trees::raster2trees() were written
  flist
  , norm_las = NULL
  , tree_sample_n = NA
  , tree_sample_prop = NA
  , which_cbh = "lowest"
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
  ##################################
  # check flist
  ##################################
    msg <- paste0(
      "If attempting to pass a list of files, the file list must:"
      , "\n   * be a vector of class character -AND-"
      , "\n   * be a directory that has final_detected_crowns* files from cloud2trees::cloud2trees() or cloud2trees::raster2trees()"
      , "\n   * -OR- be a vector of class character that includes spatial files that can be read by sf::st_read()"
    )
    if(!inherits(flist,"character")){
      stop(msg)
    }
    # check if is dir and look for cloud2trees files in the dir
    if(
      length(flist) == 1
      && dir.exists(flist)
    ){
      search_dir_final_detected_ans <- search_dir_final_detected(flist)
      crowns_flist <- search_dir_final_detected_ans$crowns_flist
      ttops_flist <- search_dir_final_detected_ans$ttops_flist
      if(is.null(crowns_flist)){stop(msg)}
    }else{
      crowns_flist <- flist
      ttops_flist <- NULL
    }
  ##################################
  # check_trees_poly
  ##################################
    # check_trees_poly will throw error if fails any checks
    check_trees_poly_ans <- crowns_flist %>% purrr::map(check_trees_poly)
  ##################################
  # sample
  ##################################
    if(
      length(ttops_flist)==length(crowns_flist)
    ){
      trees_sample <- sample_trees_flist(
        flist = ttops_flist
        , tree_sample_n = tree_sample_n
        , tree_sample_prop = tree_sample_prop
      )
    }else{
      trees_sample <- sample_trees_flist(
        flist = crowns_flist
        , tree_sample_n = tree_sample_n
        , tree_sample_prop = tree_sample_prop
      )
    }
  ##################################
  # map over trees_cbh_sf
  ##################################
  cbh_flist <-
    crowns_flist %>%
    purrr::map(function(x){
      message(paste0(
        "Attempting to extract CBH on\n   "
        , x
        , "\n   started at..."
        , Sys.time()
      ))
      ans <- trees_cbh_sf(
        trees_poly = x
        , norm_las = norm_las
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
        , trees_sample = trees_sample
        , ofile = T
      )
      return(ans)
    })
  cbh_flist <- unlist(cbh_flist)
  # return
  return(cbh_flist)
}

################################################################################################################################################
# function to sample from a list of files
# if the tree polygon data is split across multiple files
# pass this function the list of files with spatial tree crown polygons
# and it will return sample from all polygons found in the list of files
################################################################################################################################################
sample_trees_flist <- function(
  flist
  , tree_sample_n = NA
  , tree_sample_prop = NA
) {
  ##################################
  # check sample proportion
  ##################################
    check_sample_vals_ans <- check_sample_vals(
      tree_sample_n = tree_sample_n
      , tree_sample_prop = tree_sample_prop
      , def_tree_sample_n = 333
    )
    tree_sample_n <- check_sample_vals_ans$tree_sample_n
    tree_sample_prop <- check_sample_vals_ans$tree_sample_prop
  ##################################
  # check flist
  ##################################
    msg <- paste0(
      "If attempting to pass a list of files, the file list must:"
      , "\n   * be a vector of class character -AND-"
      , "\n   * be a directory that has final_detected_tree_tops* or final_detected_crowns* files from cloud2trees::cloud2trees() or cloud2trees::raster2trees()"
      , "\n   * -OR- be a vector of class character that includes spatial files that can be read by sf::st_read()"
    )
    if(!inherits(flist,"character")){
      stop(msg)
    }else{
      flist <- flist %>% normalizePath() %>% unique()
    }
    # check if is dir and look for cloud2trees files in the dir
    if(
      length(flist) == 1
      && dir.exists(flist)
    ){
      search_dir_final_detected_ans <- search_dir_final_detected(flist)
      crowns_flist <- search_dir_final_detected_ans$crowns_flist
      ttops_flist <- search_dir_final_detected_ans$ttops_flist
      if(is.null(crowns_flist) && is.null(ttops_flist)){stop(msg)}
      flist <- ifelse(is.null(ttops_flist), crowns_flist, ttops_flist)
    }
  ##################################
  # read file list
  ##################################
    trees_sf <-
      flist %>%
      normalizePath() %>%
      purrr::map(function(x){
        sf::st_read(
          dsn = x
          , quiet = T
        ) %>%
        # we only need the id
        sf::st_drop_geometry() %>%
        dplyr::mutate(
          source_filename = x
        )
      }) %>%
      dplyr::bind_rows()
  ##################################
  # check for treeID
  ##################################
    ## check for cols
    check_df_cols_all_missing(
      trees_sf
      , col_names = c("treeID")
      , all_numeric = F
    )
  ##################################
  # sample
  ##################################
    trees_sf <- trees_sf %>%
      dplyr::select(treeID,source_filename) %>%
      dplyr::group_by(treeID) %>%
      dplyr::summarise(source_filename = dplyr::first(source_filename)) %>%
      dplyr::ungroup()
    # sample
    if(
      !is.na(tree_sample_prop)
      && tree_sample_prop<1
    ){
      samp_trees <- trees_sf %>%
        dplyr::slice_sample(
          prop = tree_sample_prop
        )
    }else if(
      !is.na(tree_sample_n)
      && tree_sample_n<nrow(trees_sf)
    ){
      samp_trees <- trees_sf %>%
        dplyr::slice_sample(
          n = tree_sample_n
        )
    }else{
      samp_trees <- trees_sf
    }

  # return
    return(samp_trees)
}

################################################################################################################################################
# function to check a sf object or path of a spatial file as a character object
# for the right stuff needed to run the CBH things
################################################################################################################################################
check_trees_poly <- function(fnm, orig_fnm = "the sf object") {
  if(length(fnm) > 1){"check_trees_poly() can only do one thing at a time"}
  # geometry types from sf::st_geometry_type()
  bad_sf_types <- c(
    # "GEOMETRY","POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION" # we'll let these slide
    "POINT", "LINESTRING", "MULTIPOINT", "MULTILINESTRING"
    , "CIRCULARSTRING", "COMPOUNDCURVE", "CURVEPOLYGON", "MULTICURVE"
    , "MULTISURFACE", "CURVE", "SURFACE", "POLYHEDRALSURFACE"
    , "TIN", "TRIANGLE"
  )
  if(inherits(fnm,"character")){
    # is it a dir?
    if(dir.exists(fnm)){
      stop("check_trees_poly() can only test `sf` class objects or character objects with the path to a spatial file")
    }
    # does it even exist
    if(!file.exists(fnm)){
      stop(paste0(
        "file named \n   "
        , fnm
        , " \n   not found. does it even exist?"
      ))
    }
    # get layer name
    lyr_df <- sf::st_layers(fnm)
    # any badd types?
    has_bad <- any(toupper(lyr_df$geomtype) %in% bad_sf_types)
    if(has_bad==T){
      stop(paste0(
        "geometry type of "
        , toupper(lyr_df$geomtype[1])
        , " found in \n   "
        , fnm
        , " \n   only POLYGON type allowed for `trees_poly`"
      ))
    }
    # what if type is not read using sf::st_layers?
    test_sf <- sf::st_read(
        dsn = fnm
        , query = paste("select * from", lyr_df$name, "limit 3")
        , quiet = T
      )
  }else if(inherits(fnm,"sf")){
    test_sf <- fnm
    if(is.null(orig_fnm) || is.na(orig_fnm)){
      fnm <- "the sf object"
    }else{
      fnm <- orig_fnm
    }

  }else{
    stop("check_trees_poly() can only test `sf` class objects or character objects with the path to a spatial file")
  }
  # is it even sf?
  if(!inherits(test_sf,"sf")){
    stop(paste0(
      "this file is not even a spatial file \n   "
      , fnm
    ))
  }
  # check types
  if( !(sf::st_is(test_sf, type = c("POLYGON", "MULTIPOLYGON")) %>% all()) ){
    stop(paste0(
      "non-POLYGON found in \n   "
      , fnm
      , " \n   only POLYGON type allowed for `trees_poly`"
    ))
  }
  # check column names
  f <- test_sf %>% names() %>% dplyr::coalesce("")
  if(
    !(stringr::str_equal(f, "treeID") %>% any())
  ){
    stop(paste0(
      "`trees_poly` data must contain `treeID` column to estimate missing values."
      , "\nProvide the `treeID` as a unique identifier of individual trees."
      , "\n  ", fnm
    ))
  }else{
    # check for duplicate treeID
    if(
      nrow(test_sf) != length(unique(test_sf$treeID))
    ){
      stop(paste0(
        "Duplicates found in the treeID column. Please remove duplicates and try again."
        , "\n  ", fnm
      ))
    }
    # check that treeID is numeric or character
    if(
      !inherits(test_sf$treeID, "character")
      && !inherits(test_sf$treeID, "numeric")
    ){
      stop(paste0(
        "`trees_poly` data must contain `treeID` column of class numeric or character."
        , "\nProvide the `treeID` as a unique identifier of individual trees."
        , "\n  ", fnm
      ))
    }
  }
  # check for tree_height_m
  if(
    !(stringr::str_equal(f, "tree_height_m") %>% any())
  ){
    stop(paste0(
      "`trees_poly` data must contain `tree_height_m` column to estimate missing values."
      , "\nRename the height column if it exists and ensure it is in meters."
      , "\n  ", fnm
    ))
  }
  # success if you've made it this far
  return(TRUE)
}
################################################################################################################################################
# function to search a directory for final_detected_tree_tops* and final_detected_crowns* files
################################################################################################################################################
search_dir_final_detected <- function(dir) {
  # check for crowns
  crowns_flist <- list.files(
    normalizePath(dir)
    , pattern = "final_detected_crowns.*\\.gpkg$"
    , full.names = T
  )
  if(identical(crowns_flist, character(0))){
    crowns_flist <- NULL
  }
  # check for treetops
  ttops_flist <- list.files(
    normalizePath(dir)
    , pattern = "final_detected_tree_tops.*\\.gpkg$"
    , full.names = T
  )
  if(identical(ttops_flist, character(0))){
    ttops_flist <- NULL
  }
  return(list(
    crowns_flist = crowns_flist
    , ttops_flist = ttops_flist
  ))
}
################################################################################################################################################
# function to check the which_cbh
################################################################################################################################################
check_which_cbh <- function(which_cbh = NA) {
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
  return(which_cbh)
}
################################################################################################################################################
# function to check the tree_sample_n and tree_sample_prop
################################################################################################################################################
check_sample_vals <- function(
  tree_sample_n
  , tree_sample_prop
  , def_tree_sample_n = 333
  , def_tree_sample_prop = 0.5
) {
if(
    is.na(as.numeric(tree_sample_n)) && is.na(as.numeric(tree_sample_prop))
  ){
    tree_sample_n <- def_tree_sample_n
  }else if(
    !is.na(as.numeric(tree_sample_n)) && !is.na(as.numeric(tree_sample_prop))
  ){
    tree_sample_n <- dplyr::case_when(
      as.numeric(tree_sample_n)<=0 ~ def_tree_sample_n
      , T ~ as.numeric(tree_sample_n)
    )
    tree_sample_prop <- NA
  }else if(
    is.na(as.numeric(tree_sample_n)) && !is.na(as.numeric(tree_sample_prop))
  ){
    tree_sample_prop <- dplyr::case_when(
      as.numeric(tree_sample_prop)<=0 ~ def_tree_sample_prop
      , as.numeric(tree_sample_prop)>1 ~ 1
      , T ~ as.numeric(tree_sample_prop)
    )
    tree_sample_n <- NA
  }else if(
    !is.na(as.numeric(tree_sample_n)) && is.na(as.numeric(tree_sample_prop))
  ){
    tree_sample_n <- dplyr::case_when(
      as.numeric(tree_sample_n)<=0 ~ def_tree_sample_n
      , T ~ as.numeric(tree_sample_n)
    )
    tree_sample_prop <- NA
  }else{
    tree_sample_n <- def_tree_sample_n
    tree_sample_prop <- NA
  }
  # return
  return(list(
    tree_sample_n = tree_sample_n
    , tree_sample_prop = tree_sample_prop
  ))
}

################################################################################################################################################
## function to clip the point cloud to a polygon
## and run it through:
## leafr_for_ladderfuelsr() which is the only
## step in the process that uses the las
## and returns a data.frame
## to pass to ladderfuelsr_cbh()
################################################################################################################################################
ctg_leafr_for_ladderfuelsr <- function(
  chunk
  , poly_df
  , force_crs = F
  , voxel_grain_size_m = 1
){
  las <- lidR::readLAS(chunk)
  if(lidR::is.empty(las)){return(NULL)}
  # check for treeID in las already
  #### !!! this causes issues for if treeID in the las is not the same as in poly:
    #### !!! 1) trees don't match for checking cbh vs height
    #### !!! 2) trees don't match for estimating missing cbh values
    #### !!! how to make this work if treeID is alread in las?
  #### if(names(las@data) %>% stringr::str_equal("treeID") %>% any()){
  ####   nlas_tree <- las
  #### }else{
    # attach treeID
    nlas_tree <- polygon_attribute_to_las(
      las
      , simplify_multipolygon_crowns(poly_df) # if already simplified, does nothing
      , attribute = "treeID"
      , force_crs = force_crs
    )
  #### }
  # get the lad profile for each treeID
  safe_leafr_for_ladderfuelsr <- purrr::safely(leafr_for_ladderfuelsr)
  lad_profile <- safe_leafr_for_ladderfuelsr(
      nlas_tree
      , voxel_grain_size_m = voxel_grain_size_m
      , k = 1
      , attribute = "treeID"
      , min_pulses = 6
      , relative = F
    )
  # just get the result
  lad_profile <- lad_profile$result
  # return
  return(lad_profile)
}
################################################################################################################################################
## function to handle the lidR::catalog_apply
## output and run it through:
## ladderfuelsr_cbh() and pull out cbh values
################################################################################################################################################
output_ctg_for_ladderfuelsr_cbh <- function(
  output = NULL
  , id_class = NULL
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
      dplyr::filter(!is.na(treeID)) %>%
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

  # cast treeID in original type
  if(!inherits(lad_profile$treeID, id_class)){
    if(id_class=="character"){
      lad_profile <- lad_profile %>%
        dplyr::mutate(treeID = as.character(treeID))
    }
    if(id_class=="numeric"){
      lad_profile <- lad_profile %>%
        dplyr::mutate(treeID = as.numeric(treeID))
    }
  }

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
  # return
  return(lad_profile)
}

################################################################################################################################################
## function to check the
## return data from
## ladderfuelsr_cbh
################################################################################################################################################
clean_cbh_df <- function(cbh_df = NULL, trees_poly, lad_profile, force_cbh_lte_ht, which_cbh) {
  if(
    inherits(cbh_df, "data.frame")
    && dplyr::coalesce(nrow(cbh_df),0)>0
  ){
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

    # rid dups
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
    cbh_df <- trees_poly %>%
      dplyr::select(treeID,tree_height_m) %>%
      dplyr::inner_join(cbh_df, by = "treeID") %>%
      dplyr::filter(
        !is.na(tree_cbh_m)
        & tree_cbh_m < tree_height_m
      ) %>%
      dplyr::mutate(is_training_cbh=T) %>%
      make_spatial_predictors()
  }else if(nrow(cbh_df)>0){
    # filter cbh
    cbh_df <- trees_poly %>%
      dplyr::select(treeID) %>%
      dplyr::inner_join(cbh_df, by = "treeID") %>%
      dplyr::filter(
        !is.na(tree_cbh_m)
      ) %>%
      dplyr::mutate(is_training_cbh=T) %>%
      make_spatial_predictors()
  }

  # return
  return(cbh_df)
}

################################################################################################################################################
## function to make spatial predictor variables from POLYGON sf data with treeID and tree_height_m
# prep data for missing data model by creating predictor vars suffixed "_zzz"
################################################################################################################################################
make_spatial_predictors <- function(poly_sf) {
  df <- poly_sf %>%
    # prep data for missing data model by creating predictor vars suffixed "_zzz"
    dplyr::mutate(crown_area_zzz = sf::st_area(.) %>% as.numeric()) %>%
    sf::st_centroid() %>%
    dplyr::mutate(
      tree_x_zzz = sf::st_coordinates(.)[,1]
      , tree_y_zzz = sf::st_coordinates(.)[,2]
    ) %>%
    sf::st_drop_geometry()
  return(df)
}
