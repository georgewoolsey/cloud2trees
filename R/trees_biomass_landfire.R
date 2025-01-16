#' @title Estimate tree biomass for a tree list using LANDFIRE data
#'
#' @description
#' `trees_biomass_landfire()` uses the input tree list (e.g. as exported by [raster2trees()]) with the columns
#' `treeID`, `tree_x`, `tree_y` to attach an estimate of tree crown biomass using LANDFIRE's Forest Canopy Bulk Density (CBD) data
#' produced jointly by the U.S. Department of Agriculture and U.S. Department of the Interior.
#' If a spatial data frame of points is the input tree list, then the columns `tree_x`, `tree_y` are not required.
#'
#' LANDFIRE's Forest Canopy Bulk Density (CBD) data is attached to each tree in the tree list based on the spatial overlap with the raster data set (see references).
#' Canopy Bulk Density is mass of flammable material per unit volume of the tree crown typically expressed in units of mass per unit volume (e.g., kilograms per cubic meter).
#'
#' The simplified process for attaching tree crown biomass in kilograms to a tree is:
#'
#' * Nearest neighbor imputation is used to fill LANDFIRE data if a tree falls inside a non-forest cell in the original data
#' * The LANDFIRE estimate of CBD is distributed across the individual trees that fall in a raster cell by:
#'    1) asdf
#'    2) derp
#'    3) hey
#'
#' @param tree_list data.frame. A data frame with the columns `treeID`, `tree_x`, `tree_y`, and `tree_height_m`.
#' If an `sf` class object with POINT geometry (see [sf::st_geometry_type()]), the program will use the data "as-is" and only require the `treeID` column.
#' @param crs string. A crs string as returned from [sf::st_crs()] or the EPSG code of the x,y coordinates.
#' Defaults to the crs of the `tree_list` data if of class "sf".
#' @param study_boundary sf. The boundary of the study are to define the area of the regional model.
#' If no boundary given, regional model will be built from location of trees in the tree list.
#' @param input_landfire_dir directory where LANDFIRE Forest Canopy Bulk Density data exists. Use [get_landfire()] first.
#' @param max_search_dist_m number. Maximum search distance (m) to obtain forest type group data for trees in `tree_list` that overlap with non-forest data in the original LANDFIRE data.
#' Larger search distances will increase processing time and possibly result in memory issues.
#'
#' @references
#' * [LANDFIRE Forest Canopy Bulk Density (CBD)](https://landfire.gov/fuel/cbd)
#' U.S. Department of Agriculture and U.S. Department of the Interior.
#'
#' @return Returns a list of objects: tree_list = spatial data frame of individual trees; landfire_rast = raster of forest types in the area.
#'
#' @examples
#'  \dontrun{
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
#'  # plot the landfire_rast raster
#'  tl_type$landfire_rast %>% terra::plot()
#'  }
#' @export
#'
trees_biomass_landfire <- function(
  tree_list
  , crs = NA
  , study_boundary = NA
  , input_landfire_dir = NULL
  , max_search_dist_m = 1000
){
  # attach landfire cbd cell estimate to trees
  trees_landfire_cbd_ans <- trees_landfire_cbd(
    tree_list = tree_list
    , crs = crs
    , study_boundary = study_boundary
    , input_landfire_dir = input_landfire_dir
    , max_search_dist_m = max_search_dist_m
  )

  # the data is
  tree_list <- trees_landfire_cbd_ans$tree_list
  landfire_rast <- trees_landfire_cbd_ans$landfire_rast

  # now do the stuff in utils_biomass.R

  # return
  return(list(
    tree_list = tree_tops
    , landfire_rast = rast_ret
  ))
}
