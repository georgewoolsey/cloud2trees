#' @title function to attach polygon attribute to point cloud
#'
#' @description
#' `leafr_lad_voxels()` is a re-write of [leafR::lad.voxels()] which:
#'
#' * removes the requirement to use a file written to disk
#' * allows for the calculation of LAD by the `treeID` attribute so that don't have to pass individual tree point clouds
#' * updates to the use of the latest `lidR` functionality and removes the use of `sp` and `raster` functions
#' * updates the function to `tidy` data manipulation
#'
#' @param las an object of class LAS
#' @param poly_df an object of class sf with only POLYGON geometry (use cloud2trees::simplify_multipolygon_crowns() first)
#' @param attribute character. a data attribute in the `poly_df` that you want to spatially attach to the `las` (e.g. treeID)
#' @param force_crs logical. turn on the force same crs parameter if confident that data are in same projection.
#'
#' @return Returns an `LAS` class object
#'
#' @examples
#'  \dontrun{
#'  # polygon data
#'  f <- system.file(package = "cloud2trees","extdata","crowns_poly.gpkg")
#'  trees_poly <- sf::st_read(f)
#'  # simplify polygons
#'  trees_poly <- simplify_multipolygon_crowns(trees_poly)
#'  # point cloud data
#'  lf <- system.file(package = "cloud2trees","extdata","norm_las","RMNP_017_2018_normalize.las")
#'  las <- lidR::readLAS(lf)
#'  las@data %>% dplyr::glimpse()
#'  # polygon_attribute_to_las to attach treeID to las
#'  las <- polygon_attribute_to_las(las, trees_poly, force_crs = T, attribute = "treeID")
#'  las@data %>% dplyr::glimpse()
#'  }
#' @export
#'
polygon_attribute_to_las <- function(las, poly_df, attribute, force_crs = F){
  # check polygons
  sf_msg <- paste0(
    "`poly` data must be an object of class `sf` with only POLYGON type."
    , "\nProvide an `sf` object and see `sf::st_geometry_type()`."
  )
  if(!inherits(poly_df, "sf")){stop(sf_msg)}
  if( sf::st_is(poly_df, type = c("MULTIPOLYGON")) %>% any() ){
    stop("MULTIPOLYGON geometry detected...try cloud2trees::simplify_multipolygon_crowns() ?")
  }
  if( !(sf::st_is(poly_df, type = c("POLYGON")) %>% all()) ){stop(sf_msg)}
  # check column
  if(
    !inherits(attribute, "character") &&
    length(attribute)!=1
  ){
    stop("`attribute` must be character of length 1")
  }else(
    # check_df_cols_all_missing() in utils_biomass.r
    check_df_cols_all_missing(poly_df, col_names = attribute, all_numeric = F)
  )

  # check las data
  if (lidR::is.empty(las)) return(NULL)
  # make sure same crs
  l_epsg <- lidR::st_crs(las, parameters = T)$epsg %>%
    as.character() %>%
    dplyr::coalesce("lll")
  p_epsg <- sf::st_crs(poly_df, parameters = T) %>%
    purrr::pluck("epsg") %>%
    as.character() %>%
    dplyr::coalesce("ppp")
  if(
    force_crs &&
    !identical(lidR::st_crs(las), sf::st_crs(poly_df)) &&
    !identical(l_epsg, p_epsg)
  ){
    lidR::st_crs(las) <- sf::st_crs(poly_df)
  }else if(
    !identical(lidR::st_crs(las), sf::st_crs(poly_df)) &&
    !identical(l_epsg, p_epsg)
  ){
    stop(paste0(
      "lidR::st_crs(las) != sf::st_crs(poly_df) ensure data are same projection -or-"
      , "\nturn on the force same crs parameter if confident that data are in same projection"
    ))
  }
  # attach treeID
  nlas_tree <- lidR::merge_spatial(
    las = las
    , source = poly_df
    , attribute = attribute
  )

  return(nlas_tree)
}
#
# polygon_attribute_to_las(las, tree_crowns, force_crs = T, attribute = "treeID") %>%
#   lidR::filter_poi(!is.na(treeID)) %>%
#   lidR::plot(color = "treeID", legend = F)
