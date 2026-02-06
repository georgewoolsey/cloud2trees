#' @title Estimate HMD using tree crown polygons and normalized point cloud data
#'
#' @description
#' `trees_hmd_sf()` uses the input tree crown polygons (e.g. as exported by [raster2trees()]) with the columns
#' `treeID` and `tree_height_m` to extracting the height of the maximum crown diameter (HMD) using height normalized point cloud data (e.g. as exported by [cloud2raster()]).
#'
#' HMD is extracted directly from the height normalized point cloud by finding the height of the non-ground point farthest from the tree center (i.e. tree top).
#'
#' An early version of this process was developed by [Andrew Sanchez Meador](https://github.com/bi0m3trics).
#'
#' There are likely to be trees for which there is insufficient data in the point cloud to successfully estimate HMD. The user can elect to estimate missing HMD values which is accomplished via:
#'
#' * Attempt to extract HMD from all trees
#' * Successfully extracted HMD trees become training data used to estimate the height-HMD allometry relationship that is spatially informed using the relative tree location compared to the training data
#' * The height and location predicting HMD model built from the point cloud training data is used to predict HMD for the non-training (i.e. missing HMD) data
#'
#' @param trees_poly sf. A `sf` class object with POLYGON geometry (see [sf::st_geometry_type()]), the program will use the data "as-is" and only require the `treeID` and `tree_height_m` columns.
#' Or the path to a single spatial polygon file.
#' @param norm_las character. a directory with nomalized las files, the path of a single .laz|.las file", -or- an object of class `LAScatalog`.
#'   It is your responsibility to ensure that the point cloud is projected the same as the `trees_poly` data
#' @param tree_sample_n,tree_sample_prop numeric. Provide either `tree_sample_n`, the number of trees, or `tree_sample_prop`, the
#'   proportion of the trees to attempt to extract a HMD from the point cloud for.
#'   If neither are supplied, `tree_sample_n = 777` will be used. If both are supplied, `tree_sample_n` will be used.
#'   Increasing `tree_sample_prop` toward one (1) will increase the processing time, perhaps significantly depending on the number of trees in the `trees_poly` data.
#' @param force_same_crs logical. force the same crs between the point cloud and polygon if confident that data are in same projection.
#' data created by a `cloud2trees` pipeline (e.g. [cloud2raster()]) will always have the same projection even if not recognized by `lidR` functions
#' @param trees_sample data.frame. provide your own tree sample list such as one generated from `sample_trees_flist()` that includes the `treeID` column.
#' If provided, the tree_sample_n,tree_sample_prop will be ignored
#' @param ofile character or logical. if a character value is provided the output will be written to the disk as a csv at the location provided.
#' If set to TRUE and a file path was used as the input for `trees_poly`, then a csv file will be written to the same location with the same name prefixed with "hmd_".
#' Leave as NA to return a data.frame of the trees from tree list from `trees_poly` with HMD values added
#'
#' @noRd
#'
trees_hmd_sf <- function(
  trees_poly
  , norm_las = NULL
  , tree_sample_n = NA
  , tree_sample_prop = NA
  , force_same_crs = F
  , trees_sample = NA # answer from sample_trees_flist()
  , ofile = NA
){
  # could move to parameters
  force_hmd_lte_ht = T
  ##################################
  # check sample proportion
  ##################################
    check_sample_vals_ans <- check_sample_vals(
      tree_sample_n = tree_sample_n
      , tree_sample_prop = tree_sample_prop
      , def_tree_sample_n = 777
    )
    tree_sample_n <- check_sample_vals_ans$tree_sample_n
    tree_sample_prop <- check_sample_vals_ans$tree_sample_prop
  ##################################
  # ensure that norm las data exists
  ##################################
  nlas_ctg <- check_las_data(norm_las)
  # set the lascatalog options
  if(inherits(nlas_ctg, "LAScatalog")){
    lidR::opt_progress(nlas_ctg) <- F
    lidR::opt_filter(nlas_ctg) <- "-drop_duplicates -drop_class 2 9 18" ## class 2 = ground; 9 = water; 18 = noise
    lidR::opt_select(nlas_ctg) <- "xyz"
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
        , "\n   use trees_hmd_flist() instead?"
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
        , "max_crown_diam_height_m"
        , "is_training_hmd"
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
  # apply the ctg_calc_tree_hmd function
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
          , FUN = ctg_calc_tree_hmd
          , .options = list(automerge = TRUE)
          # ctg_calc_tree_hmd options
          , poly_df = simp_trees_poly %>% dplyr::filter(grp==x)
          , force_crs = force_same_crs
        )
      }) %>%
      unlist()
  }else{
    # apply it
    output_temp <- lidR::catalog_apply(
      ctg = nlas_ctg
      , FUN = ctg_calc_tree_hmd
      , .options = list(automerge = TRUE)
      # ctg_calc_tree_hmd options
      , poly_df = simp_trees_poly
      , force_crs = force_same_crs
    )
  }

  ##################################
  # read result from calc_tree_hmd
  ##################################
  if(
    stringr::str_ends(output_temp, ".*\\.(txt|csv)$") %>% any()
  ){
    # read the output file(s)
    hmd_df <- stringr::str_subset(output_temp, pattern = ".*\\.(txt|csv)$") %>%
      purrr::map(\(x) readr::read_delim(
        file = x, progress = F, show_col_types = F
      )) %>%
      dplyr::bind_rows() %>%
      # for trees on many tiles keep row with most points
      dplyr::filter(!is.na(treeID)) %>%
      dplyr::group_by(treeID) %>%
      dplyr::filter(calc_tree_hmd_n_pts == max(calc_tree_hmd_n_pts)) %>%
      dplyr::summarise(
        calc_tree_hmd_n_pts = max(calc_tree_hmd_n_pts, na.rm = T)
        , calc_tree_hmd_max_z = max(calc_tree_hmd_max_z, na.rm = T)
        # favor lower hmd if still duplicates
        , max_crown_diam_height_m = min(max_crown_diam_height_m, na.rm = T)
      ) %>%
      dplyr::ungroup()

    # cast treeID in original type
    if(!inherits(hmd_df$treeID, id_class)){
      if(id_class=="character"){
        hmd_df <- hmd_df %>%
          dplyr::mutate(treeID = as_character_safe(treeID))
      }
      if(id_class=="numeric"){
        hmd_df <- hmd_df %>%
          dplyr::mutate(treeID = as.numeric(treeID))
      }
    }

    # join to original data
    trees_poly <- trees_poly %>%
      dplyr::inner_join(
        hmd_df %>%
          dplyr::mutate(
            max_crown_diam_height_m = ifelse(calc_tree_hmd_n_pts<5, as.numeric(NA), max_crown_diam_height_m)
          ) %>%
          dplyr::select(-tidyselect::starts_with("calc_tree_hmd"))
        , by = "treeID"
      )
  }else{
    return(NULL)
  }

  # check force_hmd_lte_ht
  if(
    force_hmd_lte_ht==T &&
    (names(trees_poly) %>% stringr::str_equal("tree_height_m") %>% any())
  ){
    trees_poly <- trees_poly %>%
      dplyr::mutate(
        max_crown_diam_height_m = dplyr::case_when(
          is.na(max_crown_diam_height_m) ~ as.numeric(NA)
          , max_crown_diam_height_m > tree_height_m ~ as.numeric(NA)
          , T ~ max_crown_diam_height_m
        )
      )
  }

  # flag is_training_hmd
  trees_poly <- trees_poly %>%
    dplyr::mutate(
      is_training_hmd = !is.na(max_crown_diam_height_m)
    ) %>%
    dplyr::filter(is_training_hmd==T) %>%
    dplyr::select(treeID, is_training_hmd, tree_height_m, max_crown_diam_height_m) %>%
    make_spatial_predictors()

  ##################################
  # check if write it and return
  ##################################
  # anything?
  if(nrow(trees_poly)==0){return(NULL)}
  if(
    !is.na(ofile)
    && !is.null(ofile)
    && ofile==T
    && !is.null(fn_for_ofile)
    && !is.na(fn_for_ofile)
  ){
    fnm <- paste0(
      "hmd_"
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
    trees_poly %>%
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
      trees_poly %>%
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
        , "hmd_tree_list.csv"
      )
      # write it
      trees_poly %>%
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
      trees_poly %>%
        sf::st_drop_geometry() %>%
        write.csv(file = fpth, row.names = F, append = F)
      # return
      return(fpth)
    }
  }else{
    # return
    return(trees_poly)
  }
}

################################################################################################################################################
# function to apply trees_hmd_sf to a directory where cloud2trees::cloud2trees() output was written
# or to a list of crown files
# this function pipeline was developed to avoid memory issues with xxl tree lists
################################################################################################################################################
trees_hmd_flist <- function(
  # flist == a list of crown files or the directory where final_detected_tree_tops* and final_detected_crowns* files
  # from cloud2trees::cloud2trees() or cloud2trees::raster2trees() were written
  flist
  , norm_las = NULL
  , tree_sample_n = NA
  , tree_sample_prop = NA
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
  # map over trees_hmd_sf
  ##################################
  hmd_flist <-
    crowns_flist %>%
    purrr::map(function(x){
      message(paste0(
        "Attempting to extract HMD on\n   "
        , x
        , "\n   started at..."
        , Sys.time()
      ))
      ans <- trees_hmd_sf(
        trees_poly = x
        , norm_las = norm_las
        , force_same_crs = force_same_crs
        , trees_sample = trees_sample
        , ofile = T
      )
      return(ans)
    })
  hmd_flist <- unlist(hmd_flist)
  # return
  return(hmd_flist)
}

################################################################################################################################################
# let's make a function to ingest LAS class data or a data frame
  # with `x`, `y`, `z` and `treeID` columns and return the data aggregated
  # to the tree level with the HMD value
################################################################################################################################################
calc_tree_hmd <- function(las, id=NULL) {
  #######################
  # check data type
  #######################
    las_msg <- paste0(
      "`las` must contain a data frame -or- an object of class `LAS`"
      , "\nPlease update the `las` parameter."
    )
    # check
    if(inherits(las, "LAS")){
      dta <- las@data
    }else if(inherits(las,"data.frame")){
      dta <- las
    }else{stop(las_msg)}

  #######################
  # check columns
  #######################
  # overwrite treeID if id is filled
  if(!is.null(id)){
    dta <- dta %>% dplyr::mutate(treeID = id)
  }
  # names
  nms <- names(dta) %>% dplyr::coalesce("")
  # check for treeID column
  if(
    !any(stringr::str_equal(nms, "treeID"))
  ){
    stop("the `las` data does not contain the column `treeID`, ensure this column exists or set the `id` parameter")
  }
  # check for xyz
  has_xyz <- c("x", "y", "z") %>%
    purrr::map(function(x){
      stringr::str_equal(tolower(nms), x) %>%
      max() # do any columns match, T=1
    }) %>%
    unlist() %>%
    min()
  if(has_xyz==0){
    stop("the `las` data does not contain the columns `x`, `y`, and `z`, ensure columns exist")
  }

  #######################
  # find the center and farthest point
  #######################
  # classify points
  pts_temp <- dta %>%
    dplyr::rename_with(.cols = -c(treeID), .fn = tolower) %>%
    dplyr::mutate(dplyr::across(
      .cols = c(x,y,z)
      , .fns = as.numeric
    )) %>%
    dplyr::filter(!is.na(x) & !is.na(y) & !is.na(z)) %>%
    dplyr::group_by(treeID) %>% # this is key
    # first, let's arrange the points by distance from x,y center
    dplyr::mutate(
      x_mean = mean(x, na.rm = T)
      , y_mean = mean(y, na.rm = T)
      , dist_mean = sqrt((x - x_mean)^2 + (y - y_mean)^2)
    ) %>%
    # now find the point farthest from the tree center as defined by the max z point
    dplyr::mutate(
      # find the highest point in the tree as the tree "center"
      max_z = max(z, na.rm = T)
      , is_tree_center = z==max_z
      # first in case many points have max z
      , tree_center_x = dplyr::first(
        ifelse(is_tree_center, x, NA)
        , order_by = dist_mean # tie breaker is point closest to xy center
        , na_rm = T
      )
      # first in case many points have max z
      , tree_center_y = dplyr::first(
        ifelse(is_tree_center, y, NA)
        , order_by = dist_mean # tie breaker is point closest to xy center
        , na_rm = T
      )
      # Calculate the distance from the center to the point
      , dist_to_center = sqrt((x - tree_center_x)^2 + (y - tree_center_y)^2)
      , is_max_dist_to_center = dist_to_center == max(dist_to_center)
      , n_pts = dplyr::n()
    ) %>%
    dplyr::ungroup()

  # aggregate to treeID level
  df_r <- pts_temp %>%
    dplyr::filter(is_max_dist_to_center) %>%
    dplyr::group_by(treeID) %>%
    # min to be conservative for fire models
    # "conservative" = fire models will say that fire is "worse" than could be
    dplyr::summarise(
      calc_tree_hmd_max_z = dplyr::first(max_z, na_rm = T) # already aggregated, could go in group_by
      , calc_tree_hmd_n_pts = dplyr::first(n_pts, na_rm = T) # already aggregated, could go in group_by
      , max_crown_diam_height_m = min(z, na.rm = T)
    ) %>%
    dplyr::ungroup()

  # ensure that max_crown_diam_height_m<=max_z
  # return
  return(df_r)
}
################################################################################################################################################
## apply the calc_tree_hmd function to the lascatalog
################################################################################################################################################
ctg_calc_tree_hmd <- function(chunk, poly_df, force_crs = F){
  las <- lidR::readLAS(chunk)
  if (lidR::is.empty(las)) return(NULL)
  # attach treeID
  nlas_tree <- polygon_attribute_to_las(
    las
    , poly_df
    , attribute = "treeID"
    , force_crs = force_crs
  )
  # calc_tree_hmd()
  df <- calc_tree_hmd(nlas_tree)
  # return
  # return(nlas_tree)
  return(df)
}
