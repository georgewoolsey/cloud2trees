#' @title Extract LANDFIRE CBD raster cell value for a tree list based on location
#'
#' @description
#' `trees_landfire_cbd()` uses the input tree list (e.g. as exported by [raster2trees()]) with the columns
#' `treeID`, `tree_x`, `tree_y` to attach LANDFIRE's Forest Canopy Bulk Density (CBD) data estimate in kilograms per cubic meter
#' produced jointly by the U.S. Department of Agriculture and U.S. Department of the Interior.
#' If a spatial data frame of points is the input tree list, then the columns `tree_x`, `tree_y` are not required.
#'
#' LANDFIRE's Forest Canopy Bulk Density (CBD) data is attached to each tree in the tree list based on the spatial overlap with the raster data set (see references).
#' Canopy Bulk Density is mass of flammable material per unit volume of the tree crown typically expressed in units of mass per unit volume (e.g., kilograms per cubic meter).
#'
#' The simplified process for attaching the LANDFIRE CBD raster cell value to a tree is:
#'
#' * Nearest neighbor imputation is used to fill LANDFIRE data if a tree falls inside a non-forest cell in the original data
#' * The LANDFIRE raster cell value in kilograms per cubic meter is applied to a tree based on spatial overlap
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
#' * [LANDFIRE Forest Canopy Bulk Density (CBD)](https://landfire.gov/fuel/cbd)
#' U.S. Department of Agriculture and U.S. Department of the Interior.
#'
#' @return Returns a list of objects: tree_list = spatial data frame of individual trees with the column `landfire_cell_kg_per_m3` added; landfire_rast = raster of kilograms per cubic meter in the area.
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
#'  tl_lf <- trees_landfire_cbd(tree_list = tl, crs = "32613")
#'  # what?
#'  tl_lf %>% class()
#'  # a list, but what is in it?
#'  tl_lf %>% names()
#'  # what's in the trees data?
#'  tl_lf$tree_list %>% dplyr::glimpse()
#'  # plot the tree_list spatial points
#'  tl_lf$tree_list %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=landfire_cell_kg_per_m3))
#'  # plot the landfire cbd raster
#'  tl_lf$landfire_rast %>% terra::plot()
#'  # we can overlay these
#'  tl_lf$landfire_rast %>%
#'    terra::as.data.frame(xy = T) %>%
#'    ggplot2::ggplot() +
#'    ggplot2::geom_tile(ggplot2::aes(x=x,y=y,fill=kg_per_m3), color = "gray") +
#'    ggplot2::geom_sf(
#'      data = tl_lf$tree_list %>% sf::st_transform(terra::crs(tl_lf$landfire_rast))
#'      , mapping = ggplot2::aes(color=landfire_cell_kg_per_m3)
#'    ) +
#'    ggplot2::scale_color_distiller(palette = "Oranges")
#'  }
#' @export
#'
trees_landfire_cbd <- function(
  tree_list
  , crs = NA
  , study_boundary = NA
  , input_landfire_dir = NULL
  , max_search_dist_m = 1000
){
  ####################################################################
  # check external data
  ####################################################################
    # find external data
    find_ext_data_ans <- find_ext_data(
      input_landfire_dir = input_landfire_dir
    )
    # if can't find external landfire data
    if(is.null(find_ext_data_ans$landfire_dir)){
      stop(paste0(
        "LANDFIRE CBD data has not been downloaded to package contents. Use `get_landfire()` first."
        , "\nIf you supplied a value to the `input_landfire_dir` parameter check that directory for data."
      ))
    }
  ##################################
  # convert to spatial points data
  ##################################
  tree_tops <- check_spatial_points(tree_list, crs)

  # get rid of columns we'll create
    tree_tops <- tree_tops %>%
      # throw in hey_xxxxxxxxxx to test it works if we include non-existant columns
      dplyr::select( -dplyr::any_of(c(
        "hey_xxxxxxxxxx"
        , "landfire_cell_kg_per_m3"
      )))

  ##################################
  # load landfire data. see get_landfire()
  ##################################
    # get the landfire data
    landfire <- terra::rast(
        file.path(find_ext_data_ans$landfire_dir, "lc23_cbd_240.tif")
      )

  ####################################################################
  # crop the raster and extract values at point locations
  ####################################################################
    # call crop_raster_match_points() defined in utils_rast_points.R
    crop_raster_match_points_ans <- crop_raster_match_points(
      points = tree_tops
      , rast = landfire
      , study_boundary = study_boundary
      , max_search_dist_m = max_search_dist_m
    )

    ###############################
    # define return data
    ###############################
    point_values_from_rast <- crop_raster_match_points_ans$point_values
    rast_ret <- crop_raster_match_points_ans$rast

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for landfire we just need to ensure that we get >0, non-NA values
    # let's check with the values
    na_trees <- crop_raster_match_points_ans$point_values %>%
      dplyr::filter(
        is.na(as.numeric(raster_value))
        | tolower(raster_value) %in% c("fill-nodata", "non-forested")
        | stringr::str_detect(tolower(raster_value), "nodata")
        | stringr::str_detect(tolower(raster_value), "non")
        | as.numeric(raster_value)<=0
      ) %>%
      nrow()

  ####################################################################
  # IF we even need to fill NA values
  ####################################################################
  if(na_trees>0){
    # mark all non-desired cells as "NA" in the cropped raster
    ### non-desired cells change depending on the raster data used !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # reclass_landfire_rast() defined in utils_rast_points.R
     reclass_rast <- reclass_landfire_rast(
       rast = crop_raster_match_points_ans$rast
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
          landfire_cell_kg_per_m3 = point_values_from_rast$raster_value %>% as.numeric()
        )
    }else{
      tree_tops <- tree_tops %>%
        dplyr::mutate(
          landfire_cell_kg_per_m3 = as.numeric(NA)
        )
      message(paste0(
        "Unable to determine LANDFIRE CBD for this tree list and study boundary (if provided)."
        , "\nTry expanding the study boundary area or increasing the max_search_dist_m parameter"
        , "\nand ensure that your tree data is in the continental US."
      ))
    }

  # rename the return raster
    terra::set.names(rast_ret, value = "kg_per_m3")

  # return
  return(list(
    tree_list = tree_tops
    , landfire_rast = rast_ret
  ))
}
