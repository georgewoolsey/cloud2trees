#' @title Simplify MULTIPOLYGON to POLYGON geometry in an `sf` class object
#'
#' @param trees_poly data.frame. A data frame of `sf` class with POLYGON,MULTIPOLYGON geometry (see [sf::st_geometry_type()]) and the column `treeID`
#'
#' @description
#' Function to simplify MULTIPOLYGON geometry to POLYGON geometry in an `sf` class object by selecting the largest segment of the MULTIPOLYGON
#'
#' @return A `sf` class object data frame
#'
#' @examples
#'  \dontrun{
#'  f <- paste0(system.file(package = "cloud2trees"),"/extdata/crowns_poly.gpkg")
#'  crowns <- sf::st_read(f, quiet = T)
#'  crowns %>% sf::st_geometry_type() %>% table()
#'  crowns_simp <- simplify_multipolygon_crowns(crowns)
#'  crowns_simp %>% sf::st_geometry_type() %>% table()
#'  }
#' @export
#'
simplify_multipolygon_crowns <- function(trees_poly) {
  if(!inherits(trees_poly, "sf")){
    stop("must pass an sf object to the `trees_poly` parameter")
  }

  # check if not polygon
  if( min(sf::st_is(trees_poly, type = c("POLYGON", "MULTIPOLYGON"))) == 0 ){
    warning(paste0(
      "data passed to `trees_poly` is not polygon or multipolygon data"
      , "\n see sf::st_geometry_type...returning original data"
    ))
    return(trees_poly)
  }

  # return data as-is if no multipolygon
  if( max(sf::st_is(trees_poly, type = c("MULTIPOLYGON"))) == 0 ){
    return(trees_poly)
  }

  # check if has a treeID
  if(
    (names(trees_poly) %>% stringr::str_detect("treeID") %>% max())==0
  ){
    stop(paste0(
      "`trees_poly` data must contain `treeID` column."
      , "\nProvide the `treeID` as a unique identifier of individual trees."
    ))
  }

  # simplify the multipolygons
  df <-
    # start with only polygons
    trees_poly %>%
    dplyr::filter(sf::st_geometry_type(.)=="POLYGON") %>%
    # union on cleaned multipolygons
    dplyr::bind_rows(
      trees_poly %>%
        dplyr::filter(sf::st_geometry_type(.)=="MULTIPOLYGON") %>%
        sf::st_cast(to = "POLYGON", do_split = T, warn = F) %>%
        dplyr::mutate(axxx = sf::st_area(.)) %>% # axxx is so we don't overwrite a column
        dplyr::group_by(treeID) %>%
        dplyr::filter(axxx == max(axxx)) %>% # keep the biggest crown polygon by treeID
        dplyr::ungroup() %>%
        dplyr::select(-axxx)
    ) %>%
    dplyr::ungroup()

  # return
  return(df)
}
