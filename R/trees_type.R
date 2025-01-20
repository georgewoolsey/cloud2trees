#' @title Estimate forest type for a tree list based on location
#'
#' @description
#' `trees_type()` uses the input tree list (e.g. as exported by [raster2trees()]) with the columns
#' `treeID`, `tree_x`, `tree_y` to attach species information using USDA Forest Inventory and Analysis (FIA) codes.
#' If a spatial data frame of points is the input tree list, then the columns `tree_x`, `tree_y` are not required.
#'
#' FIA Forest Type Group Code is attached to each tree in the tree list based on the spatial overlap with the Forest Type Groups of the Continental United States dataset [Wilson 2023](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d).
#'
#' The simplified process for attaching forest type group to a tree is:
#'
#' * Forest type group 30-m raster (Wilson 2023) was aggregated to 90-m to make the data more accessible over the entire continental US
#' * Nearest neighbor imputation is used to fill forest type data if a tree falls inside a no-forest cell in the original data
#' * The FIA forest type group is applied to a tree based on spatial overlap
#'
#' @param tree_list data.frame. A data frame with the columns `treeID`, `tree_x`, `tree_y`, and `tree_height_m`.
#' If an `sf` class object with POINT geometry (see [sf::st_geometry_type()]), the program will use the data "as-is" and only require the `treeID` column.
#' @param crs string. A crs string as returned from [sf::st_crs()] or the EPSG code of the x,y coordinates.
#' Defaults to the crs of the `tree_list` data if of class "sf".
#' @param study_boundary sf. The boundary of the study are to define the area of the regional model.
#' If no boundary given, regional model will be built from location of trees in the tree list.
#' @param input_foresttype_dir directory where Forest Type Groups data exists. Use [get_foresttype()] first.
#' @param max_search_dist_m number. Maximum search distance (m) to obtain forest type group data for trees in `tree_list` that overlap with non-forest data in the original Wilson (2023) data.
#' Larger search distances will increase processing time and possibly result in memory issues.
#'
#' @references
#' * [Forest Type Groups of the Continental United States](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d)
#' Wilson, B.T. (2023). Forest Type Groups of the Continental United States.
#'
#' @return Returns a list of objects: tree_list = spatial data frame of individual trees; foresttype_rast = raster of forest types in the area.
#'
#' @examples
#'  \dontrun{
#'  library(tidyverse)
#'  # example tree list
#'  tl <- dplyr::tibble(
#'      treeID = c(1:21)
#'      , tree_x = rnorm(n=21, mean = 458064, sd = 11)
#'      , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
#'    )
#'  # call the function
#'  tl_type <- trees_type(tree_list = tl, crs = "32613")
#'  # what?
#'  tl_type %>% class()
#'  # a list, but what is in it?
#'  tl_type %>% names()
#'  # plot the tree_list spatial points
#'  tl_type$tree_list %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=forest_type_group))
#'  # plot the foresttype_rast raster
#'  tl_type$foresttype_rast %>% terra::plot()
#'  }
#' @export
#'
trees_type <- function(
  tree_list
  , crs = NA
  , study_boundary = NA
  , input_foresttype_dir = NULL
  , max_search_dist_m = 1000
){
  ####################################################################
  # check external data
  ####################################################################
    # find external data
    find_ext_data_ans <- find_ext_data(
      input_foresttype_dir = input_foresttype_dir
    )
    # if can't find external foresttype data
    if(is.null(find_ext_data_ans$foresttype_dir)){
      stop(paste0(
        "Forest Type Group data has not been downloaded to package contents. Use `get_foresttype()` first."
        , "\nIf you supplied a value to the `input_foresttype_dir` parameter check that directory for data."
      ))
    }
  ##################################
  # convert to spatial points data
  ##################################
  f <- tree_list %>% names() %>% dplyr::coalesce("") # leaving in case anything below looks for this
  tree_tops <- check_spatial_points(tree_list, crs)

  # get rid of columns we'll create
    tree_tops <- tree_tops %>%
      # throw in hey_xxxxxxxxxx to test it works if we include non-existant columns
      dplyr::select( -dplyr::any_of(c(
        "hey_xxxxxxxxxx"
        , "forest_type_group_code"
        , "forest_type_group"
        , "hardwood_softwood"
      )))

  ##################################
  # load foresttype data. see get_foresttype()
  ##################################
    # get the foresttype data
    foresttype <- terra::rast(
        file.path(find_ext_data_ans$foresttype_dir, "foresttype.tif")
      )
    # get the lookup
    foresttype_lookup <- readr::read_csv(
        file.path(find_ext_data_ans$foresttype_dir, "foresttype_lookup.csv")
        , progress = F
        , show_col_types = F
      ) %>%
      dplyr::mutate(dplyr::across(
        dplyr::everything()
        , as.character
      ))

  ####################################################################
  # crop the raster and extract values at point locations
  ####################################################################
    # call crop_raster_match_points() defined in utils_rast_points.R
    crop_raster_match_points_ans <- crop_raster_match_points(
      points = tree_tops
      , rast = foresttype
      , study_boundary = study_boundary
      , max_search_dist_m = max_search_dist_m
    )

    ###############################
    # define return data
    ###############################
    point_values_from_rast <- crop_raster_match_points_ans$point_values
    rast_ret <- crop_raster_match_points_ans$rast

    # # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for forest type we need to check with the lookup table
    # # let's check with the lookup table
    na_trees <- crop_raster_match_points_ans$point_values %>%
      dplyr::mutate(
        forest_type_code = raster_value %>% as.character()
      ) %>%
      dplyr::left_join(foresttype_lookup, by = "forest_type_code") %>%
      dplyr::filter(is.na(forest_type_group_code)) %>%
      nrow()

  ####################################################################
  # IF we even need to fill NA values
  ####################################################################
  if(na_trees>0){
    # mark all non-desired cells as "NA" in the cropped raster
    ### non-desired cells change depending on the raster data used !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # reclass_foresttype_rast() defined in utils_rast_points.R
     reclass_rast <- reclass_foresttype_rast(
       rast = crop_raster_match_points_ans$rast
       , lookup = foresttype_lookup
      )

    # aggregate (if needed) and fill missing raster values
    # agg_fill_rast_match_points() defined in utils_rast_points.R
    agg_fill_rast_match_points_ans <- agg_fill_rast_match_points(
      rast = reclass_rast
      , bbox = crop_raster_match_points_ans$bbox
      , points = tree_tops
    )

    ###############################
    # update the return data
    ###############################
    if(!is.null(agg_fill_rast_match_points_ans$point_values)){
      point_values_from_rast <- agg_fill_rast_match_points_ans$point_values
    }
    if(!is.null(agg_fill_rast_match_points_ans$rast)){
      rast_ret <- agg_fill_rast_match_points_ans$rast
    }

  }

  #######################################################
  # prep final data
  #######################################################
    if(length(point_values_from_rast$raster_value)==nrow(tree_tops)){
      # now let's join it back with our data and check it
      tree_tops <- tree_tops %>%
        dplyr::mutate(
          forest_type_code = point_values_from_rast$raster_value %>% as.character()
        ) %>%
        dplyr::left_join(
          foresttype_lookup %>%
            dplyr::select(
              forest_type_code
              , forest_type_group_code
              , forest_type_group
              , hardwood_softwood
            )
          , by = "forest_type_code"
        ) %>%
        dplyr::select(-forest_type_code)
    }else{
      tree_tops <- tree_tops %>%
        dplyr::mutate(
          forest_type_group_code = as.character(NA)
          , forest_type_group = as.character(NA)
          , hardwood_softwood = as.character(NA)
        )
      message(paste0(
        "Unable to determine forest type for this tree list and study boundary (if provided)."
        , "\nTry expanding the study boundary area or increasing the max_search_dist_m parameter"
        , "\nand ensure that your tree data is in the continental US."
      ))
    }

  # return
  return(list(
    tree_list = tree_tops
    , foresttype_rast = rast_ret
  ))
}
