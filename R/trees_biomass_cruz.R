#' @title Estimate tree crown biomass for a tree list based on Cruz et al. (2003)
#'
#' @description
#' `trees_biomass_cruz()` uses the input tree list (e.g. as exported by [raster2trees()]) with the columns
#' `treeID`, `tree_x`, `tree_y` to attempt to attach tree crown biomass in kilogram estimates based on
#' the Cruz et al. (2003) equations (see references) and the FIA forest type group.
#'
#' If a spatial data frame of points is the input tree list, then the columns `tree_x`, `tree_y` are not required.
#' If the FIA forest type group data named `forest_type_group_code`is not in the input tree list, then the function calls
#' [trees_type()] to attempt to attach the FIA forest type group.
#' Other required columns include:
#'
#' * `crown_area_m2`, `tree_height_m` (e.g. as exported by [raster2trees()])
#' * `tree_cbh_m` (e.g. as exported by [trees_cbh()])
#' * and one of `dbh_cm`, `dbh_m`, or  `basal_area_m2` (e.g. as exported by [trees_dbh()])
#'
#' The Cruz et al. (2003) study developed models to predict canopy fuel stratum at the
#' stand level for four coniferous forest types common in the western US:
#' Douglas-fir, ponderosa pine, lodgepole pine, and mixed conifer.
#' Models for other forests types are currently lacking which limits the scope of this methodology.
#' If the tree list has trees that are in a FIA forest type group not represented in the list above,
#' then the return data will be blank
#'
#' Canopy Bulk Density is mass of flammable material per unit volume of the tree
#' crown typically expressed in units of mass per unit volume (e.g., kilograms per cubic meter).
#'
#' The process for estimating tree crown biomass in kilograms is:
#'
#' * Nearest neighbor imputation is used to fill the FIA forest type data if a tree falls inside a non-forest cell in the original data
#' * The LANDFIRE estimate of CBD is distributed across the individual trees that fall in a raster cell by:
#'    1) At the stand level (i.e. raster cell), aggregate the tree level data within the stand to obtain:
#'      - `mean_crown_length_m = mean(crown_length_m)`, where tree `crown_length_m = tree_height_m - tree_cbh_m`
#'      - `sum_crown_volume_m3 = sum(crown_volume_m3)`, where tree `crown_volume_m3 = (4/3) * pi * ((crown_length_m/2)) * ((crown_dia_m/2)^2)`
#'    2) At the stand level (i.e. raster cell), determine the area of the stand that overlaps (`overlap_area_m2`) with the AOI defined as the `study_boundary` parameter (see below) or the bounding box of all the trees
#'    3) At the stand level (i.e. raster cell), use the Cruz equations (Table 4; see reference) to estimate of CBD in kilograms per cubic meter (`cruz_stand_kg_per_m3`)
#'    4) At the stand level (i.e. raster cell), get canopy fuel loading (CFL) in kilograms per square meter (`kg_per_m2 = mean_crown_length_m * cruz_stand_kg_per_m3`)
#'    5) At the stand level (i.e. raster cell), get the stand biomass in kilograms (`biomass_kg = kg_per_m2 * overlap_area_m2`)
#'    6) At the stand level (i.e. raster cell), the single tree CBD in kilograms per cubic meter will be a constant (`cruz_tree_kg_per_m3 = biomass_kg / sum_crown_volume_m3`)
#'    7) Attach the the single tree CBD in kilograms per cubic meter to the tree level based on raster cell spatial overlap
#'    8) Calculate individual tree crown mass in kilograms as `cruz_crown_biomass_kg = cruz_tree_kg_per_m3 * crown_volume_m3`
#'
#' @param tree_list data.frame. A data frame with the columns `treeID`, `tree_x`, and `tree_y`.
#' If an `sf` class object with POINT geometry (see [sf::st_geometry_type()]), the program will use the data "as-is" and only require the `treeID` column.
#' Other required columns include:
#' * `crown_area_m2`, `tree_height_m` (e.g. as exported by [raster2trees()])
#' * `tree_cbh_m` (e.g. as exported by [trees_cbh()])
#' * and one of `dbh_cm`, `dbh_m`, or  `basal_area_m2` (e.g. as exported by [trees_dbh()])
#' @param crs string. A crs string as returned from [sf::st_crs()] or the EPSG code of the x,y coordinates.
#' Defaults to the crs of the `tree_list` data if of class "sf".
#' @param study_boundary sf. The boundary of the study area to define the area of interest which may extend beyond the space with trees.
#' If no boundary given, the AOI will be built from location of trees in the tree list.
#' @param input_foresttype_dir directory where Forest Type Groups data exists. Use [get_foresttype()] first.
#' @param max_crown_kg_per_m3 numeric. the maximum CBD of the tree crown in kilograms per cubic meter.
#' Values above this limit will be set at the median value for the area using only stands that have CBD values lower than this limit.
#' The default value of 2 kilograms per cubic meter was based on [Mell et al. (2009)](https://doi.org/10.1016/j.combustflame.2009.06.015)
#' who found the dry bulk density of the tree crown was 2.6 kilograms per cubed meter
#' using Douglas-fir trees grown on Christmas tree farms.
#' Set this parameter to a large value (e.g. 1e10) or NULL to avoid limiting tree crown CBD.
#'
#' @references
#' * [Forest Type Groups of the Continental United States](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d)
#' Wilson, B.T. (2023). Forest Type Groups of the Continental United States.
#' * [doi:10.1071/WF02024](https://scholar.google.com/scholar?cluster=316241498622221569&oi=gsb&hl=en&as_sdt=0,5)
#' Cruz, M.G, M.E. Alexander, and R.H. Wakimoto. 2003. Assessing canopy fuel stratum characteristics in crown fire prone fuel types of western North America. Int. J. Wildland Fire. 12(1):39-50.
#' * [https://doi.org/10.1016/j.combustflame.2009.06.015](https://doi.org/10.1016/j.combustflame.2009.06.015)
#' Mell, W., Maranghides, A., McDermott, R., & Manzello, S. L. (2009). Numerical simulation and experiments of burning douglas fir trees. Combustion and Flame, 156(10), 2023-2041.
#'
#'
#' @return Returns a list of objects: tree_list = spatial data frame of individual trees
#' ; stand_cell_data = data frame of stands/cells in same projection as the FIA forest type group raster data
#'
#' See code in examples.
#'
#' @examples
#'  \dontrun{
#' library(tidyverse)
#' library(sf)
#' my_n <- 111
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
#' tl_cruz <- trees_biomass_cruz(tree_list = tl, crs = "32613")
#' # what is in it?
#' tl_cruz %>% names()
#' # look at the trees
#' tl_cruz$tree_list %>% dplyr::glimpse()
#' # tree FIA forest type groups
#' tl_cruz$tree_list %>%
#'   sf::st_drop_geometry() %>%
#'   dplyr::count(forest_type_group_code, forest_type_group)
#' # look at the stand
#' tl_cruz$stand_cell_data %>% dplyr::filter(!is.na(trees)) %>% dplyr::glimpse()
#' # get the projection for the stand cell data
#' epsg_code <- tl_cruz$stand_cell_data$rast_epsg_code[1] %>% as.numeric()
#'  # plot the stand cell data with trees overlaid
#'  tl_cruz$stand_cell_data %>%
#'    ggplot2::ggplot() +
#'    ggplot2::geom_tile(ggplot2::aes(x=x,y=y,fill = cruz_stand_kg_per_m3), color = "gray44") +
#'    ggplot2::geom_text(ggplot2::aes(x=x,y=y,label = trees), color = "white") +
#'    ggplot2::geom_sf(
#'      data = tl_cruz$tree_list %>% sf::st_transform(crs = epsg_code)
#'      , ggplot2::aes(color = cruz_crown_biomass_kg)
#'    ) +
#'    ggplot2::labs(fill="stand kg/m3", color = "tree crown kg", caption = "# trees shown in cell") +
#'    ggplot2::scale_fill_viridis_c(option = "rocket", na.value = "gray", direction = -1) +
#'    ggplot2::scale_color_viridis_c(option = "viridis", na.value = "gray22", begin = 0.6) +
#'    ggplot2::theme_void()
#'  }
#' @export
#'
trees_biomass_cruz <- function(
  tree_list
  , crs = NA
  , study_boundary = NA
  , input_foresttype_dir = NULL
  , max_crown_kg_per_m3 = 2
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
      clean_biomass_cols(method = "cruz")

  ##################################
  # ensure that we have forest_type_group_code
  # this is required for cruz estimation
  ##################################
  nms <- tree_tops %>% names() %>% dplyr::coalesce("")
  if(
    !(stringr::str_equal(nms, "forest_type_group_code") %>% any())
  ){
    warning(paste0(
      "`tree_list` data must contain `forest_type_group_code` column."
      , "\nAttempting to determine forest type group using trees_type()..."
    ))
    safe_trees_type <- purrr::safely(trees_type)
    trees_type_ans <- safe_trees_type(
      tree_list = tree_tops
      , study_boundary = study_boundary
      , input_foresttype_dir = find_ext_data_ans$foresttype_dir
    )
    # check error
    if(is.null(trees_type_ans$error)){
      # keep going
      tree_tops <- trees_type_ans$result$tree_list
    }else{
      stop(paste0(
        "Error: failed to get forest type group using trees_type(). Message:"
        , "\n"
        , trees_type_ans$error
      ))
    }
  }
  ##################################
  # check forest_type_group_code
  # only certain types are available for cruz
  ##################################
  has_cruz <- has_cruz_forest_type_group_code(tree_tops)
  if(!has_cruz){
    warning(paste0(
        "None of the forest types present match with the Cruz equations available..."
        , "\n returning original data with forest type group if it was attached via trees_type()"
      ))
    return(tree_tops)
  }

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
        , as_character_safe
      ))

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
      rast = foresttype
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
      , cbd_method = "cruz"
      , max_crown_kg_per_m3 = max_crown_kg_per_m3
    )

  # return
  return(list(
    tree_list = distribute_stand_fuel_load_ans$tree_list
    , stand_cell_data = distribute_stand_fuel_load_ans$cell_df
  ))
}
