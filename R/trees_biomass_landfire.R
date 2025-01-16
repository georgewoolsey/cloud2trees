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
#'    1) get cfl in kg/m2 == kg_per_m2 = mean_crown_length_m * landfire_stand_kg_per_m3
#'    2) get stand biomass in kg at the stand level == biomass_kg = kg_per_m2 * overlap_area_m2
#'    3) single tree CBD in kg/m3 will be constant by stand/cell == landfire_tree_kg_per_m3 = biomass_kg / sum_crown_volume_m3
#'    4) at the tree level landfire_tree_biomass_kg = landfire_tree_kg_per_m3*crown_volume_m3
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
#' library(tidyverse)
#' library(sf)
#' my_n <- 122
#' # fake tree list
#' tl <- dplyr::tibble(
#'     treeID = c(1:my_n)
#'     , tree_x = rnorm(n=my_n, mean = 458064, sd = 33)
#'     , tree_y = rnorm(n=my_n, mean = 4450074, sd = 33)
#'     , tree_height_m = rnorm(n=my_n, mean = 10, sd = 7)
#'   ) %>%
#'   dplyr::mutate(
#'     tree_height_m = ifelse(tree_height_m<1.37, 1.37, tree_height_m) # above DBH
#'     , crown_area_m2 = 0.47+(0.49*tree_height_m)
#'     , tree_cbh_m = 0.72+(0.53*tree_height_m)
#'     , dbh_cm = -2.3+(2.14*tree_height_m)
#'   )
#' # how does our fake tree list look?
#' tl %>% dplyr::glimpse()
#' tl %>% ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(x = dbh_cm, y = tree_height_m))
#' # call the function
#' tl_landfire <- trees_biomass_landfire(tree_list = tl, crs = "32613")
#' # what is in it?
#' tl_landfire %>% names()
#' # look at the trees
#' tl_landfire$tree_list %>% dplyr::glimpse()
#' # look at the stand
#' tl_landfire$stand_cell_data %>% dplyr::filter(!is.na(trees)) %>% dplyr::glimpse()
#' # the stand CBD
#'  tl_landfire$stand_cell_data %>%
#'    dplyr::filter(trees>0) %>%
#'    sf::st_drop_geometry() %>%
#'    dplyr::count(landfire_stand_kg_per_m3)
#' # get the projection for the stand cell data
#' epsg_code <- tl_landfire$stand_cell_data$rast_epsg_code[1] %>% as.numeric()
#'  # plot the stand cell data with trees overlaid
#'  tl_landfire$stand_cell_data %>%
#'    ggplot2::ggplot() +
#'    ggplot2::geom_tile(ggplot2::aes(x=x,y=y,fill = landfire_stand_kg_per_m3), color = "gray44") +
#'    ggplot2::geom_text(ggplot2::aes(x=x,y=y,label = trees), color = "white") +
#'    ggplot2::geom_sf(
#'      data = tl_landfire$tree_list %>% sf::st_transform(crs = epsg_code)
#'      , ggplot2::aes(color = landfire_tree_biomass_kg)
#'    ) +
#'    ggplot2::labs(fill="stand kg/m3", color = "tree kg", caption = "# trees shown in cell") +
#'    ggplot2::scale_fill_viridis_c(option = "rocket", na.value = "gray", direction = -1) +
#'    ggplot2::scale_color_viridis_c(option = "viridis", na.value = "gray22", begin = 0.6) +
#'    ggplot2::theme_void()
#'  }
#' @export
#'
trees_biomass_landfire <- function(
  tree_list
  , crs = NA
  , study_boundary = NA
  , input_landfire_dir = NULL
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
  # check columns needed for biomass estimation
  check_df_cols_all_missing(
      tree_tops
      , col_names = c("crown_area_m2", "tree_height_m", "tree_cbh_m")
      , all_numeric = T
    )
  # check for DBH or BA
  safe_check_df_cols_all_missing <- purrr::safely(check_df_cols_all_missing)
  ba_chk <- c("dbh_cm","dbh_m","basal_area_m2") %>%
    purrr::map(\(x)
      safe_check_df_cols_all_missing(col_names = x, df = tree_tops, all_numeric = T)
    ) %>%
    purrr::transpose() %>%
    purrr::pluck("result") %>%
    purrr::flatten() %>%
    unlist()
  if(!any(ba_chk)){
    stop(paste0(
        "the data does not contain the columns `basal_area_m2`, `dbh_cm`, or `dbh_m`"
        , "\n .... at least one of these columns must be present and have data"
      ))
  }

  # get rid of columns we'll create
    tree_tops <- tree_tops %>%
      # throw in hey_xxxxxxxxxx to test it works if we include non-existent columns
      dplyr::select( -dplyr::any_of(c(
        "hey_xxxxxxxxxx"
        , "landfire_stand_id"
        , "crown_dia_m"
        , "crown_length_m"
        , "crown_volume_m3"
        , "landfire_tree_kg_per_m3"
        , "landfire_stand_kg_per_m3"
        , "landfire_tree_biomass_kg"
      )))

  ##################################
  # ensure that we have landfire_cell_kg_per_m3
  # this is required for landfire estimation
  ##################################
  nms <- tree_tops %>% names() %>% dplyr::coalesce("")
  if(
    !(stringr::str_equal(nms, "landfire_cell_kg_per_m3") %>% any())
  ){
    warning(paste0(
      "`tree_list` data must contain `landfire_cell_kg_per_m3` column."
      , "\nAttempting to determine CBD using trees_landfire_cbd()..."
    ))
    safe_trees_landfire_cbd <- purrr::safely(trees_landfire_cbd)
    trees_landfire_cbd_ans <- safe_trees_landfire_cbd(
      tree_list = tree_tops
      , study_boundary = study_boundary
      , input_landfire_dir = find_ext_data_ans$landfire_dir
    )
    # check error
    if(is.null(trees_landfire_cbd_ans$error)){
      # keep going
      tree_tops <- trees_landfire_cbd_ans$result$tree_list
    }else{
      stop(paste0(
        "Error: failed to get LANDFIRE CBD using trees_landfire_cbd(). Message:"
        , "\n"
        , trees_landfire_cbd_ans$error
      ))
    }
  }
  ##################################
  # check landfire_cell_kg_per_m3
  ##################################
  lf_chk <- safe_check_df_cols_all_missing(col_names = "landfire_cell_kg_per_m3", df = tree_tops, all_numeric = T)
  if(!is.null(lf_chk$error)){
    stop(paste0(
      "Error: failed to get LANDFIRE CBD using trees_landfire_cbd()."
      , "\n Is this area in the continental United States? Try to expand search area in trees_landfire_cbd()?"
    ))
  }

  ##################################
  # load landfire data. see get_landfire()
  ##################################
    # get the landfire data
    landfire <- terra::rast(
        file.path(find_ext_data_ans$landfire_dir, "lc23_cbd_240.tif")
      )

  ##################################
  # create extent if empty
  ##################################
    # set study boundary to tree extent if missing
    if(
      !(
        c(
          inherits(study_boundary, "sf")
          , inherits(study_boundary, "sfc")
          , inherits(study_boundary, "SpatVector")
        ) %>% any()
      )
    ){
      study_boundary <- tree_tops %>%
        sf::st_bbox() %>%
        sf::st_as_sfc()
    }

  ################################################################
  # calc_rast_cell_trees
    # function to aggregate tree list to the raster cell level
    # and join to the raster cell overlap data generated via
    # calc_rast_cell_overlap() within the function
  ################################################################
    calc_rast_cell_trees_ans <- calc_rast_cell_trees(
      rast = landfire
      , tree_list = tree_tops
      , poly_extent = study_boundary
      , calc_tree_level_cols = T
    )

  ################################################################
  # distribute_stand_fuel_load
    # use our calculate stand-level CBD for cruz (lf already has CBD)
    # the stand level CBD in kilograms per cubed meter is distributed to trees
  ################################################################
    distribute_stand_fuel_load_ans <- distribute_stand_fuel_load(
      cell_df = calc_rast_cell_trees_ans$cell_df
      , tree_list = calc_rast_cell_trees_ans$tree_list
      , cbd_method = "landfire"
    )

  # return
  return(list(
    tree_list = distribute_stand_fuel_load_ans$tree_list
    , stand_cell_data = distribute_stand_fuel_load_ans$cell_df
  ))
}
