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
#' @param tree_sample_n,tree_sample_prop numeric. Provide either `tree_sample_n`, the number of trees, or `tree_sample_prop`, the
#'   proportion of the trees to attempt to extract a HMD from the point cloud for.
#'   If neither are supplied, `tree_sample_n = 777` will be used. If both are supplied, `tree_sample_n` will be used.
#'   Increasing `tree_sample_prop` toward one (1) will increase the processing time, perhaps significantly depending on the number of trees in the `trees_poly` data.
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
  , tree_sample_n = NA
  , tree_sample_prop = NA
  , estimate_missing_hmd = F
  , force_same_crs = F
){
  # could move to parameters
  force_hmd_lte_ht = T

  ##############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  xxx

  # ensure that there are enough data to estimate
  n_hmd <- trees_poly %>%
    sf::st_drop_geometry() %>%
    dplyr::filter(is_training_hmd==T) %>%
    nrow()

  #######################################################
  # estimate_missing_hmd
  #######################################################
  n_hmd <- dplyr::coalesce(n_hmd,0)
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
      # predictors and response to pass to randomForest functions
      predictors <- training_df %>% dplyr::select(-c(treeID,max_crown_diam_height_m))
      response <- training_df$max_crown_diam_height_m

      # implements steps to mitigate very long run-times when tuning random forests models
      optimal_mtry <- rf_tune_subsample(
        predictors = predictors
        , response = response
      )

      ### Run a randomForest model to predict HMD using various crown predictors
      # quiet this
      quiet_rf <- purrr::quietly(randomForest::randomForest)
      # run it
      hmd_mod <- quiet_rf(
        y = response
        , x = predictors
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

  }else if(n_hmd==0){
    message(paste0(
      "No HMD values extracted"
    ))
    return(trees_poly)
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
