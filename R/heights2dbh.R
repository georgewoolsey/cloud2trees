#' @title Use a CHM raster to detect individual trees
#'
#' @description
#' `heights2dbh()` is an all-in-one function to process a CHM raster
#' and return a spatial data frame of tree crown polygons and points.
#' The order of operations is:
#'
#' * Perform individual tree detection using [lidR::locate_trees()] with the [lidR::lmf()] algorithm
#' * Delineate tree crowns using [ForestTools::mcws()]
#'
#' Note, this function does not estimate DBH for the detected trees and only returns tree location, crown area, and height information.
#' To estimate tree DBH from the detected tree heights see [heights2dbh()].
#'
#' @param chm_rast raster. A  raster from `terra` or `stars`representing a canopy height model
#' @param outfolder string. The path of a folder to write the crown vector data to
#' @param ws numeric or function. Length or diameter of the moving window used to detect the local
#' maxima in the units of the input data (usually meters). If it is numeric a fixed window size is used.
#' If it is a function, the function determines the size of the window at any given location on the canopy.
#' By default function takes the height of a given pixel as its only argument and return the
#' desired size of the search window when centered on that pixel.
#' @param min_height numeric. Set the minimum height (m) for individual tree detection
#' @param min_crown_area numeric. Set the minimum crown area (m2) for individual tree detection
#' @param tempdir string. Directory to write intermediate files. Intermediate files are only created for large rasters too big to fit in memory.
#'
#' @references
#' https://r-lidar.github.io/lidRbook/itd.html
#'
#' @return Returns a spatial data frame of individual tree crown vectors detected using the CHM.
#' The tree top point coordinates are located in the `tree_x` and `tree_y` columns.
#' The process also writes two `.gpkg` files to the `outfolder` directory: `chm_detected_crowns.gpkg` and `chm_detected_tree_tops.gpkg`
#'
#' @examples
#'  \dontrun{
#'  o <- "../data"
#'  i <- "../data/lasdata"
#'  r <- cloud2trees::cloud2raster(output_dir = o, input_las_dir = i)
#'  r %>% names()
#'  r$dtm_rast %>% terra::plot()
#'  r$chm_rast %>% terra::plot()
#'  r$create_project_structure_ans %>% dplyr::glimpse()
#'  r$chunk_las_catalog_ans$process_data %>% dplyr::glimpse()
#'  r$chunk_las_catalog_ans$grid_subset_switch
#'  r$chunk_las_catalog_ans$las_ctg@data %>% dplyr::glimpse()
#'  r$normalize_flist
#'  }
#' @export
#'
heights2dbh <- function(x) {
  return(x)
}
