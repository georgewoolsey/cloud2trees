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
#' @param trees_poly sf. A `sf` class object with POLYGON geometry (see [sf::st_geometry_type()]), the program will use the data "as-is" and only require the `treeID` and `tree_height_m` columns.
#' @param norm_las character. a directory with nomalized las files, the path of a single .laz|.las file", -or- an object of class `LAScatalog`.
#'   It is your responsibility to ensure that the point cloud is projected the same as the `trees_poly` data
#' @param estimate_missing_hmd logical. it is not likely that HMD will be extracted successfully from every tree (especially in low density clouds).
#'   Should the missing HMD values be estimated using the tree height and location information based on trees for which HMD is successfully extracted?
#' @param force_same_crs logical. force the same crs between the point cloud and polygon if confident that data are in same projection.
#' data created by a `cloud2trees` pipeline (e.g. [cloud2raster()]) will always have the same projection even if not recognized by `lidR` functions
#'
#' @references
#' An early version of this process was developed by [Andrew Sanchez Meador](https://github.com/bi0m3trics).
#'
#' @return Returns a spatial data frame of individual trees with the added columns: `max_crown_diam_height_m`, `is_training_hmd`
#'
#' @examples
#'  \dontrun{
#'  # example tree crown polygons
#'  f <- paste0(system.file(package = "cloud2trees"),"/extdata/crowns_poly.gpkg")
#'  crowns <- sf::st_read(f, quiet = T)
#'  # example normalized las files are in this directory
#'  norm_d <- paste0(system.file(package = "cloud2trees"),"/extdata/norm_las")
#'  # now run the trees_hmd()
#'  trees_hmd_ans <- trees_hmd(
#'     trees_poly = crowns
#'     , norm_las = norm_d
#'     , force_same_crs = T
#'     , estimate_missing_hmd = T)
#'  # what?
#'  trees_hmd_ans %>% dplyr::glimpse()
#'  # spatial polygons
#'  trees_hmd_ans %>% ggplot2::ggplot() +
#'     ggplot2::geom_sf(ggplot2::aes(fill=max_crown_diam_height_m))
#'  # relationship between height and hmd
#'  trees_hmd_ans %>%
#'     ggplot2::ggplot(
#'       ggplot2::aes(
#'         x = tree_height_m, y = max_crown_diam_height_m, color=is_training_hmd
#'       )
#'     ) +
#'     ggplot2::geom_point()
#'  }
#' @export
#'
trees_hmd <- function(
  trees_poly
  , norm_las = NULL
  , estimate_missing_hmd = F
  , force_same_crs = F
){
  # could move to parameters
  force_hmd_lte_ht = T
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
  # ensure spatial polygon data
  ##################################
  sf_msg <- paste0(
      "`trees_poly` data must be an object of class `sf` with only POLYGON type."
      , "\nProvide an `sf` object and see `sf::st_geometry_type()`."
    )
  if(!inherits(trees_poly, "sf")){stop(sf_msg)}
  if( !(sf::st_is(trees_poly, type = c("POLYGON", "MULTIPOLYGON")) %>% all()) ){stop(sf_msg)}

  ##################################
  # ensure that treeID data exists
  ##################################
  f <- trees_poly %>% names() %>% dplyr::coalesce("")
  if(
    !(stringr::str_equal(f, "treeID") %>% any())
  ){
    stop(paste0(
      "`trees_poly` data must contain `treeID` column to estimate HMD."
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
      "`trees_poly` data must contain `tree_height_m` column to estimate HMD."
      , "\nRename the height column if it exists and ensure it is in meters."
    ))
  }

  # get rid of columns we'll create
    trees_poly <- trees_poly %>%
      # throw in hey_xxxxxxxxxx to test it works if we include non-existant columns
      dplyr::select( -dplyr::any_of(c(
        "hey_xxxxxxxxxx"
        , "max_crown_diam_height_m"
        , "is_training_hmd"
      )))

  ##################################
  # apply the ctg_calc_tree_hmd function
  ##################################
  # simplify the polygons so that lidR::merge_spatial can be used
  simp_trees_poly <- simplify_multipolygon_crowns(trees_poly)

  # check if we need to split for massive tree crown data
  if(nrow(simp_trees_poly)>500e3){
    # for data with so many crowns I've encountered the error:
    ### !!! Error in getGlobalsAndPackages(expr, envir = envir, tweak = tweakExpression, :
      ### !!! The total size of the xx globals exported for future expression ...
      ### !!! is xxx GiB.. This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize').

    # apply it
    output_temp <- lidR::catalog_apply(
      ctg = nlas_ctg
      , FUN = ctg_calc_tree_hmd
      , .options = list(automerge = TRUE)
      # ctg_calc_tree_hmd options
      , poly_df = simp_trees_poly
      , force_crs = force_same_crs
    )
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
    id_class <- class(trees_poly$treeID)[1]
    if(!inherits(hmd_df$treeID, id_class)){
      if(id_class=="character"){
        lad_profile <- lad_profile %>%
          dplyr::mutate(treeID = as.character(treeID))
      }
      if(id_class=="numeric"){
        lad_profile <- lad_profile %>%
          dplyr::mutate(treeID = as.numeric(treeID))
      }
    }

    # join to original data
    trees_poly <- trees_poly %>%
      dplyr::left_join(
        hmd_df %>%
          dplyr::mutate(
            max_crown_diam_height_m = ifelse(calc_tree_hmd_n_pts<3, as.numeric(NA), max_crown_diam_height_m)
          ) %>%
          dplyr::select(-tidyselect::starts_with("calc_tree_hmd"))
        , by = "treeID"
      )
  }else{
    # blank the hmd column
    trees_poly <- trees_poly %>%
      dplyr::mutate(
        max_crown_diam_height_m = as.numeric(NA)
      )
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
    )

  # ensure that there are enough data to estimate
  n_hmd <- trees_poly %>%
    sf::st_drop_geometry() %>%
    dplyr::filter(is_training_hmd==T) %>%
    nrow()

  #######################################################
  # estimate_missing_hmd
  #######################################################
  if(
    estimate_missing_hmd==T
    && n_hmd > 10
    && (names(trees_poly) %>% stringr::str_equal("tree_height_m") %>% any())
  ){
    # add x,y to data
    mod_df <- trees_poly %>%
      dplyr::select(treeID, is_training_hmd, tree_height_m, max_crown_diam_height_m) %>%
      dplyr::mutate(crown_area_zzz = sf::st_area(.) %>% as.numeric()) %>%
      sf::st_centroid() %>%
      dplyr::mutate(
        tree_xxx = sf::st_coordinates(.)[,1]
        , tree_yyy = sf::st_coordinates(.)[,2]
        , tree_height_m = as.numeric(tree_height_m)
        , max_crown_diam_height_m = as.numeric(max_crown_diam_height_m)
      ) %>%
      sf::st_drop_geometry()
    # training versus predict data
    training_df <- mod_df %>% dplyr::filter(is_training_hmd==T) %>% dplyr::select(-is_training_hmd)
    predict_df <- mod_df %>% dplyr::filter(is_training_hmd==F) %>% dplyr::select(-is_training_hmd)

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
        y = training_df$max_crown_diam_height_m
        , x = training_df %>% dplyr::select(-c(treeID,max_crown_diam_height_m))
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
          training_df %>% dplyr::select(-c(treeID,max_crown_diam_height_m))
        )
      )

      ### Run a randomForest model to predict HMD using various crown predictors
      # quiet this
      quiet_rf <- purrr::quietly(randomForest::randomForest)
      # run it
      hmd_mod <- quiet_rf(
        y = training_df$max_crown_diam_height_m
        , x = training_df %>% dplyr::select(-c(treeID,max_crown_diam_height_m))
        , mtry = optimal_mtry
        , na.action = na.omit
      )
      # just get the result
      hmd_mod <- hmd_mod$result

    # # model
    # hmd_mod <- stats::lm(
    #   formula = max_crown_diam_height_m ~ tree_xxx + tree_yyy + tree_xxx:tree_yyy + tree_height_m + crown_area_zzz
    #   , data = training_df
    # )

    # predict missing
    predicted_hmd_temp <- predict(
        hmd_mod
        , predict_df %>% dplyr::select(-c(treeID,max_crown_diam_height_m))
      ) %>%
      dplyr::as_tibble() %>%
      dplyr::pull(1)

    ## combine predicted data with training data for full data set
    trees_poly <- trees_poly %>%
      # join with predicted data estimates
      dplyr::left_join(
        predict_df %>%
          dplyr::mutate(
            predicted_hmd = predicted_hmd_temp
          ) %>%
          dplyr::select(treeID, predicted_hmd)
        , by = dplyr::join_by("treeID")
      ) %>%
      # clean up data
      dplyr::mutate(
        max_crown_diam_height_m = dplyr::coalesce(max_crown_diam_height_m, predicted_hmd)
      ) %>%
      dplyr::select(-predicted_hmd)

  }else if(estimate_missing_cbh==T){
    if(max(grepl("tree_height_m", f))==0){
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
  }

  ## prevent the hmd from being > the tree height
    if(force_hmd_lte_ht==T && max(grepl("tree_height_m", f))==1){
      # find the 95th percentile of height-cbh ratio
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

#####################################################
#####################################################
# intermediate functions
#####################################################
#####################################################
################
# let's make a function to ingest LAS class data or a data frame
  # with `x`, `y`, `z` and `treeID` columns and return the data aggregated
  # to the tree level with the HMD value
################
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
################
## apply the function to the lascatalog
################
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
