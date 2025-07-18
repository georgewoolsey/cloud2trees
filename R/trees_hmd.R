#' @title Estimate HMD using tree crown polygons and normalized point cloud data
#'
#' @description
#' `trees_hmd()` uses the input tree crown polygons (e.g. as exported by [raster2trees()]) with the columns
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
#' @param trees_poly must be one of the following that has required attributes `treeID` and `tree_height_m`:
#'   * `sf` class object with POLYGON geometry (see [sf::st_geometry_type()]). Recommended for smaller tree lists (e.g. <100k) that can fit in memory.
#'   * character vector with the path to a single or multiple spatial files that can be read by [sf::st_read()] and have POLYGON geometry. Recommended for large tree lists (e.g. 100k+) that might cause memory issues.
#'   * character with the path to a directory that has "final_detected_crowns*" files from [cloud2trees()] or [raster2trees()]. Recommended for large tree lists (e.g. 100k+) that might cause memory issues.
#'
#' @param norm_las character. a directory with nomalized las files, the path of a single .laz|.las file", -or- an object of class `LAScatalog`.
#'   It is your responsibility to ensure that the point cloud is projected the same as the `trees_poly` data
#' @param tree_sample_n,tree_sample_prop numeric. Provide either `tree_sample_n`, the number of trees, or `tree_sample_prop`, the
#'   proportion of the trees to attempt to extract a HMD from the point cloud for.
#'   If neither are supplied, `tree_sample_n = 777` will be used. If both are supplied, `tree_sample_n` will be used.
#'   Increasing `tree_sample_prop` toward one (1) will increase the processing time, perhaps significantly depending on the number of trees in the `trees_poly` data.
#' @param estimate_missing_hmd logical. it is not likely that HMD will be extracted successfully from every tree (especially in low density clouds).
#'   Should the missing HMD values be estimated using the tree height and location information based on trees for which HMD is successfully extracted?
#' @param force_same_crs logical. force the same crs between the point cloud and polygon if confident that data are in same projection.
#' data created by a `cloud2trees` pipeline (e.g. [cloud2raster()]) will always have the same projection even if not recognized by `lidR` functions
#' @param outfolder string. The path of a folder to write the model data to.
#'   Note, in the actual missing value estimation many RF models are estimated and model averaging is used.
#'   However, only the first estimated model is saved in this export which does not fully represent the process used to fill in missing values.
#'
#' @references
#' An early version of this process was developed by [Andrew Sanchez Meador](https://github.com/bi0m3trics).
#'
#' @return Returns a spatial data frame of individual trees with the added columns: `max_crown_diam_height_m`, `is_training_hmd`
#'
#' @examples
#'  \dontrun{
#'   library(tidyverse)
#'   library(sf)
#'   # example tree crown polygons
#'   f <- paste0(system.file(package = "cloud2trees"),"/extdata/crowns_poly.gpkg")
#'   crowns <- sf::st_read(f, quiet = T)
#'   # example normalized las files are in this directory
#'   norm_d <- paste0(system.file(package = "cloud2trees"),"/extdata/norm_las")
#'   # now run the trees_hmd()
#'   trees_hmd_ans <- trees_hmd(
#'      trees_poly = crowns
#'      , norm_las = norm_d
#'      , force_same_crs = T
#'      , estimate_missing_hmd = T)
#'   # what?
#'   trees_hmd_ans %>% dplyr::glimpse()
#'   # spatial polygons
#'   trees_hmd_ans %>% ggplot2::ggplot() +
#'      ggplot2::geom_sf(ggplot2::aes(fill=max_crown_diam_height_m))
#'   # relationship between height and hmd
#'   trees_hmd_ans %>%
#'      ggplot2::ggplot(
#'        ggplot2::aes(
#'          x = tree_height_m, y = max_crown_diam_height_m, color=is_training_hmd
#'        )
#'      ) +
#'      ggplot2::geom_point()
#'   #### try a file list
#'   #### Recommended for large tree lists (e.g. 100k+) that might cause memory issues.
#'   # we'll split the crowns
#'   # as is done automatically for tree lists >250k by raster2trees() and cloud2trees()
#'   crowns <- crowns %>%
#'     dplyr::mutate(
#'       # makes 2 groups of data
#'       grp = ceiling(dplyr::row_number()/(dplyr::n()/2))
#'     )
#'   # make file names
#'   my_dir <- tempdir()
#'   fnm_1 <- file.path(my_dir, "crowns1.gpkg")
#'   fnm_2 <- file.path(my_dir, "crowns2.gpkg")
#'   fnm_1
#'   # write the data
#'   sf::st_write(crowns %>% dplyr::filter(grp==1), dsn = fnm_1, append = F) # grp 1
#'   sf::st_write(crowns %>% dplyr::filter(grp==2), dsn = fnm_2, append = F) # grp 2
#'   # try trees_cbh with our file list
#'   flist <- c(fnm_1,fnm_2)
#'   # now run the trees_hmd()
#'   trees_hmd_ans2 <- trees_hmd(
#'      trees_poly = flist
#'      , norm_las = norm_d
#'      , force_same_crs = T
#'      , estimate_missing_hmd = T)
#'   # tabulate training data
#'   trees_hmd_ans %>%
#'     sf::st_drop_geometry() %>%
#'     dplyr::count(is_training_hmd)
#'  }
#' @export
#'
trees_hmd <- function(
  trees_poly
  , norm_las = NULL
  , tree_sample_n = NA
  , tree_sample_prop = NA
  , estimate_missing_hmd = TRUE
  , force_same_crs = F
  , outfolder = tempdir()
){
  # could move to parameters
  force_hmd_lte_ht <- T
  estimate_missing_max_n_training <- 25000
  #############################################
  # estimate hmd based on what's in trees_poly
  #############################################
  crowns_flist <- NULL # setting up for filling
  if(inherits(trees_poly, "sf")){
    #############################################
    # if trees_poly is sf
    #############################################
    hmd_df <- trees_hmd_sf(
      trees_poly = trees_poly
      , norm_las = norm_las
      , tree_sample_n = tree_sample_n
      , tree_sample_prop = tree_sample_prop
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
      }else if(length(crowns_flist)==1){ # if only one file then trees_hmd_sf() so that only read file once
        hmd_df <- trees_hmd_sf(
          trees_poly = crowns_flist
          , norm_las = norm_las
          , tree_sample_n = tree_sample_n
          , tree_sample_prop = tree_sample_prop
          , force_same_crs = force_same_crs
        )
      }else{
        hmd_df <- trees_hmd_flist(
          flist = trees_poly # just searches the directory again so that tree points are used for sample
          , norm_las = norm_las
          , tree_sample_n = tree_sample_n
          , tree_sample_prop = tree_sample_prop
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
      hmd_df <- trees_hmd_sf(
        trees_poly = trees_poly
        , norm_las = norm_las
        , tree_sample_n = tree_sample_n
        , tree_sample_prop = tree_sample_prop
        , force_same_crs = force_same_crs
      )
    }else{
      #############################################
      # if trees_poly is a list of filenames
      #############################################
      crowns_flist <- trees_poly %>% normalizePath() %>% unique()
      hmd_df <- trees_hmd_flist(
        flist = crowns_flist
        , norm_las = norm_las
        , tree_sample_n = tree_sample_n
        , tree_sample_prop = tree_sample_prop
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
  # read hmd_df if it was a file list to get
  # trees where hmd extraction was successful
  #############################################
  if(
    inherits(hmd_df, "character")
    && (stringr::str_ends(hmd_df, ".*\\.csv$") %>% any())
  ){
    # read the output file(s)
    hmd_df <- stringr::str_subset(hmd_df, pattern = ".*\\.csv$") %>%
      readr::read_csv(progress = F, show_col_types = F) %>%
      dplyr::mutate(dplyr::across(
        .cols = c(tree_height_m, max_crown_diam_height_m, tidyselect::ends_with("_zzz"))
        , .fns = as.numeric
      )) %>%
      dplyr::filter(is_training_hmd==T)
  }else if(inherits(hmd_df, "data.frame")){
    hmd_df <- hmd_df %>%
      dplyr::mutate(dplyr::across(
        .cols = c(tree_height_m, max_crown_diam_height_m, tidyselect::ends_with("_zzz"))
        , .fns = as.numeric
      )) %>%
      dplyr::filter(is_training_hmd==T)
  }else{
    stop("error extracting HMD")
  }
  # ensure that there are enough data to estimate
  n_hmd <- hmd_df %>%
    sf::st_drop_geometry() %>%
    dplyr::filter(is_training_hmd==T) %>%
    nrow()

  #######################################################
  # build hmd model using training data
  # happens before we load trees_poly if inherits(crowns_flist, "character")
  # so that memory isn't already stuffed
  #######################################################
  n_hmd <- dplyr::coalesce(n_hmd,0)
  if(
    estimate_missing_hmd==T
    && n_hmd > 10
    && (names(hmd_df) %>% stringr::str_equal("tree_height_m") %>% any())
  ){
    # go for at least 50% of training data sampled
    # it's just that the individual model fits will be smaller
    ntimes_temp <- ((nrow(hmd_df)*0.5)/estimate_missing_max_n_training) %>%
      ceiling() %>%
      max(3) %>%
      min(50)
    # estimate
    hmd_mod <- rf_subsample_and_model_n_times(
      predictors = hmd_df %>% dplyr::select(-c(treeID,max_crown_diam_height_m,is_training_hmd))
      , response = hmd_df$max_crown_diam_height_m
      , mod_n_subsample = estimate_missing_max_n_training
      , mod_n_times = ntimes_temp
    )
    # write just the first model
    if(
      !is.null(hmd_mod)
      && length(hmd_mod)>0
    ){
      saveRDS(hmd_mod[[1]], file = file.path(normalizePath(outfolder), "hmd_height_model_estimates.rds"))
    }
  }else{
    hmd_mod <- NULL
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
        , "max_crown_diam_height_m"
        , "is_training_hmd"
      )))
  }else if(
    !is.null(crowns_flist)
    && inherits(crowns_flist, "character")
  ){
    # if we've made it this far, the polygon data has already gone through the checks in trees_hmd_sf()
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
          , "max_crown_diam_height_m"
          , "is_training_hmd"
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

  #######################################################
  # fill missing hmd values
  #######################################################
  if(
    estimate_missing_hmd==T
    && !is.null(hmd_mod)
    && length(hmd_mod)>0
    && (names(trees_poly) %>% stringr::str_equal("tree_height_m") %>% any())
  ){
    # model the missing values
    predict_df <- trees_poly %>%
      dplyr::anti_join(hmd_df %>% dplyr::select(treeID), by = "treeID") %>%
      dplyr::select(treeID, tree_height_m) %>%
      make_spatial_predictors()
    # get predicted values
    predicted_temp <- rf_model_avg_predictions(mod_list = hmd_mod, predict_df = predict_df)
    # attach back to data
    predict_df <- predict_df %>%
      dplyr::mutate(predicted_zzz = predicted_temp$predicted)

    ## combine predicted data with training data for full data set
    trees_poly <- trees_poly %>%
      # join with training data
      dplyr::left_join(
        hmd_df %>%
          dplyr::select(treeID, max_crown_diam_height_m, is_training_hmd)
        , by = "treeID"
      ) %>%
      dplyr::mutate(is_training_hmd = dplyr::coalesce(is_training_hmd, F)) %>%
      # join with predicted data estimates
      dplyr::left_join(
        predict_df %>%
          dplyr::select(treeID, predicted_zzz)
        , by = dplyr::join_by("treeID")
      ) %>%
      # clean up data
      dplyr::mutate(
        max_crown_diam_height_m = dplyr::coalesce(max_crown_diam_height_m, predicted_zzz)
      ) %>%
      dplyr::select(-predicted_zzz)

  }else if(n_hmd==0){
    message(paste0(
      "No HMD values extracted"
    ))
    return(
      trees_poly %>%
        dplyr::mutate(
          max_crown_diam_height_m = as.numeric(NA)
          , is_training_hmd = as.logical(NA)
        )
    )
  }else if(estimate_missing_hmd==T){
    if(!(names(trees_poly) %>% stringr::str_equal("tree_height_m") %>% any()) ){
      message(paste0(
        "`trees_poly` data must contain `tree_height_m` column to estimate HMD."
        , "\nSetting `estimate_missing_hmd=TRUE` requires this data."
        , "\nReturning HMD values extracted from cloud only."
      ))
    }else{
      message(paste0(
        "Insufficient data available to estimate missing HMD values."
        , "\nReturning HMD values extracted from cloud only."
      ))
    }
    ## combine predicted data with training data for full data set
    trees_poly <- trees_poly %>%
      dplyr::left_join(
        hmd_df %>%
          dplyr::select(treeID, max_crown_diam_height_m, is_training_hmd)
        , by = "treeID"
      ) %>%
      dplyr::mutate(is_training_hmd = dplyr::coalesce(is_training_hmd, F))
  }else{
    ## combine predicted data with training data for full data set
    trees_poly <- trees_poly %>%
      dplyr::left_join(
        hmd_df %>%
          dplyr::select(treeID, max_crown_diam_height_m, is_training_hmd)
        , by = "treeID"
      ) %>%
      dplyr::mutate(is_training_hmd = dplyr::coalesce(is_training_hmd, F))
  }

  ## prevent the hmd from being > the tree height
    if(
      force_hmd_lte_ht==T &&
      (names(trees_poly) %>% stringr::str_equal("tree_height_m") %>% any())
    ){
      # find the 95th percentile of height-hmd ratio
      max_ratio <- trees_poly %>%
        dplyr::filter(
          is_training_hmd==T
          & max_crown_diam_height_m < tree_height_m
        ) %>%
        dplyr::mutate(ratio = max_crown_diam_height_m/tree_height_m) %>%
        dplyr::pull(ratio) %>%
        stats::quantile(probs = 0.95)
      # update values
      trees_poly <- trees_poly %>%
        dplyr::mutate(
          # update training data where max_crown_diam_height_m > tree_height_m
          is_training_hmd = dplyr::case_when(
            is_training_hmd==T & max_crown_diam_height_m > tree_height_m ~ FALSE
            , T ~ is_training_hmd
          )
          # update max_crown_diam_height_m
          , max_crown_diam_height_m = dplyr::case_when(
            is_training_hmd==F & max_crown_diam_height_m/tree_height_m > max_ratio ~ max_ratio*tree_height_m
            , T ~ max_crown_diam_height_m
          )
        )
    }

  # return
  return(trees_poly)
}
