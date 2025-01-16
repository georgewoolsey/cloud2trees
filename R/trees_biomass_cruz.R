#' @title Estimate tree biomass for a tree list based on Cruz et al. (2003)
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
#' stand level for four coniferous forest types common in the western US: Douglas-fir, ponderosa pine, lodgepole pine, and mixed conifer.
#' Models for other forests types are currently lacking which limits the scope of this methodology.
#' If the tree list has trees that fall are in a FIA forest type group not represented in the list above, then the return data will be blank
#'
#' Canopy Bulk Density is mass of flammable material per unit volume of the tree crown typically expressed in units of mass per unit volume (e.g., kilograms per cubic meter).
#'
#' The process for distributing the stand-level to a tree is:
#'
#' * do it
#' * keep doing it
#' * hey
#'
#' @param tree_list data.frame. A data frame with the columns `treeID`, `tree_x`, `tree_y`, and `tree_height_m`.
#' If an `sf` class object with POINT geometry (see [sf::st_geometry_type()]), the program will use the data "as-is" and only require the `treeID` column.
#' @param crs string. A crs string as returned from [sf::st_crs()] or the EPSG code of the x,y coordinates.
#' Defaults to the crs of the `tree_list` data if of class "sf".
#' @param study_boundary sf. The boundary of the study are to define the area of the regional model.
#' If no boundary given, regional model will be built from location of trees in the tree list.
#' @param input_foresttype_dir directory where Forest Type Groups data exists. Use [get_foresttype()] first.
#'
#' @references
#' * [Forest Type Groups of the Continental United States](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d)
#' Wilson, B.T. (2023). Forest Type Groups of the Continental United States.
#' * [doi:10.1071/WF02024](https://scholar.google.com/scholar?cluster=316241498622221569&oi=gsb&hl=en&as_sdt=0,5)
#' Cruz, M.G, M.E. Alexander, and R.H. Wakimoto. 2003. Assessing canopy fuel stratum characteristics in crown fire prone fuel types of western North America. Int. J. Wildland Fire. 12(1):39-50.
#'
#' @return Returns a list of objects: tree_list = spatial data frame of individual trees; stand_cell_data = data frame of stands/cells in same projection as the foresttype raster data
#' See code in examples.
#'
#' @examples
#'  \dontrun{
#' library(tidyverse)
#' library(sf)
#' my_n <- 55
#' # fake tree list
#' tl <- dplyr::tibble(
#'     treeID = c(1:my_n)
#'     , tree_x = rnorm(n=my_n, mean = 458064, sd = 11)
#'     , tree_y = rnorm(n=my_n, mean = 4450074, sd = 11)
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
#' # look at the stand
#' tl_cruz$stand_cell_data %>% dplyr::filter(!is.na(trees)) %>% dplyr::glimpse()
#' # get the projection for the stand cell data
#' epsg_code <- tl_cruz$stand_cell_data$rast_epsg_code[1] %>% as.numeric()
#' # plot the stand cell data with trees overlaid
#' tl_cruz$stand_cell_data %>% dplyr::filter(!is.na(trees)) %>%
#'   ggplot2::ggplot() +
#'   ggplot2::geom_tile(ggplot2::aes(x=x,y=y,fill = cruz_stand_kg_per_m3)) +
#'   ggplot2::geom_text(ggplot2::aes(x=x,y=y,label = trees), color = "white") +
#'   ggplot2::geom_sf(
#'     data = tl_cruz$tree_list %>% sf::st_transform(crs = epsg_code)
#'     , ggplot2::aes(color = cruz_tree_biomass_kg)
#'   ) +
#'   ggplot2::labs(fill="stand kg/m3", color = "tree kg", caption = "# trees shown in cell") +
#'   ggplot2::scale_fill_viridis_c(option = "rocket", na.value = NA, direction = -1) +
#'   ggplot2::scale_color_viridis_c(option = "viridis", na.value = "gray22", begin = 0.6) +
#'   ggplot2::theme_void()
#'  }
#' @export
#'
trees_biomass_cruz <- function(
  tree_list
  , crs = NA
  , study_boundary = NA
  , input_foresttype_dir = NULL
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

  # get rid of columns we'll create
    tree_tops <- tree_tops %>%
      # throw in hey_xxxxxxxxxx to test it works if we include non-existent columns
      dplyr::select( -dplyr::any_of(c(
        "hey_xxxxxxxxxx"
        , "cruz_stand_id"
        , "crown_dia_m"
        , "crown_length_m"
        , "crown_volume_m3"
        , "cruz_tree_kg_per_m3"
        , "cruz_stand_kg_per_m3"
        , "cruz_tree_biomass_kg"
      )))

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
        , as.character
      ))

  ##################################
  # create extent if empty
  ##################################
    # set study boundary to tree extent if missing
    study_boundary <- NA
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

    # mapview::mapview(
    #   study_boundary %>% sf::st_as_sf()
    #   , color = "red", lwd = 3, alpha.regions = 0
    #   , layer.name = "boundary", label = FALSE, legend = FALSE, popup = FALSE
    # ) +
    #   mapview::mapview(tree_list, zcol = "tree_height_m")
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

    # calc_rast_cell_trees_ans$cell_df %>% dplyr::filter(!is.na(trees)) %>% dplyr::glimpse()
    # #
    # calc_rast_cell_trees_ans$cell_df %>% dplyr::filter(!is.na(trees)) %>%
    #   ggplot() +
    #   geom_tile(aes(x=x,y=y,fill = trees_per_ha)) +
    #   geom_text(aes(x=x,y=y,label = trees), color = "white") +
    #   geom_sf(
    #     data = calc_rast_cell_trees_ans$tree_list %>% sf::st_transform(terra::crs(foresttype))
    #     , mapping = aes(color = tree_height_m)
    #   ) +
    #   labs(fill="TPH", color = "tree ht. (m)", caption = "# trees shown in cell") +
    #   scale_fill_viridis_c(na.value = NA) +
    #   scale_color_viridis_c(option = "rocket", na.value = "black", begin = 0.55, direction = -1) +
    #   theme_void()

  ################################################################
  # distribute_stand_fuel_load
    # use our `get_cruz_stand_kg_per_m3()` function to calculate
    # the stand level CBH in kilograms per cubed meter
    # and distribute this across the tree list
  ################################################################
    distribute_stand_fuel_load_ans <- distribute_stand_fuel_load(
      cell_df = calc_rast_cell_trees_ans$cell_df
      , tree_list = calc_rast_cell_trees_ans$tree_list
      , cbd_method = "cruz"
    )

  # return
  return(list(
    tree_list = distribute_stand_fuel_load_ans$tree_list
    , stand_cell_data = distribute_stand_fuel_load_ans$cell_df
  ))
}
