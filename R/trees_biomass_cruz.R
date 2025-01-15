#' @title Estimate tree biomass for a tree list
#'
#' @description
#' `trees_biomass()` uses the input tree list (e.g. as exported by [raster2trees()]) with the columns
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
#'
#' @references
#' * [Forest Type Groups of the Continental United States](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d)
#' Wilson, B.T. (2023). Forest Type Groups of the Continental United States.
#'
#' @return Returns a list of objects: tree_list = spatial data frame of individual trees; foresttype_rast = raster of forest types in the area.
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
#'  tl_type <- trees_biomass(tree_list = tl, crs = "32613")
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
trees_biomass <- function(
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
        , "cell"
        , "crown_dia_m"
        , "crown_length_m"
        , "crown_volume_m3"
        , "tree_kg_per_m3"
        , "cruz_biomass_kg"
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
    # use our `get_cruz_stand_kg_per_m3()` function to calculate
    # the stand level CBH in kilograms per cubed meter
    # and distribute this across the tree list
  ################################################################
    distribute_stand_fuel_load_ans <- distribute_stand_fuel_load(
      cell_df = calc_rast_cell_trees_ans$cell_df
      , tree_list = calc_rast_cell_trees_ans$tree_list
    )

  # return
  return(
    distribute_stand_fuel_load_ans$tree_list
  )
}
