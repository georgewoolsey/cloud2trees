#' @title Estimate tree biomass (or crown biomass) for a tree list
#'
#' @description
#' `trees_biomass()` streamlines the process for estimating individual tree biomass in kilograms, or the component biomass of the tree crown in kilograms.
#' Users can select one or all of the following methods available in the package for estimating biomass:
#'
#' * Tree crown biomass in kilograms:
#'    - "landfire" - based on LANDFIRE's Forest Canopy Bulk Density (CBD) data ([trees_biomass_landfire()])
#'    - "cruz" - based on Cruz et al. (2003) canopy fuel stratum equations ([trees_biomass_cruz()])
#' * Tree total above ground biomass in kilograms:
#'    - coming soon
#'
#' If multiple methods are selected (e.g. `method = c("cruz","landfire")`), then the program will compile the
#' biomass estimates and return one tree list.
#'
#' @inheritParams trees_biomass_landfire
#' @param input_foresttype_dir directory where Forest Type Groups data exists. Use [get_foresttype()] first.
#' @param method character. one (e.g. `"landfire"`) or multiple (e.g. `c("cruz","landfire")`) of the following biomass estimation methods:
#' * Tree crown biomass in kilograms:
#'    - "landfire" - based on LANDFIRE's Forest Canopy Bulk Density (CBD) data ([trees_biomass_landfire()])
#'    - "cruz" - based on Cruz et al. (2003) canopy fuel stratum equations ([trees_biomass_cruz()])
#' * Tree total above ground biomass in kilograms:
#'    - coming soon
#'
#' @references
#' See references in:
#' * [trees_biomass_landfire()]
#' * [trees_biomass_cruz()]
#'
#' @return Returns a list of objects: tree_list = spatial data frame of individual trees
#' ; stand_cell_data_landfire = data frame of stands/cells in same projection as the LANDFIRE raster data
#' ; stand_cell_data_cruz = data frame of stands/cells in same projection as the FIA forest type group raster data
#'
#' See code in examples.
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
#'      , ggplot2::aes(color = landfire_crown_biomass_kg)
#'    ) +
#'    ggplot2::labs(fill="stand kg/m3", color = "tree kg", caption = "# trees shown in cell") +
#'    ggplot2::scale_fill_viridis_c(option = "rocket", na.value = "gray", direction = -1) +
#'    ggplot2::scale_color_viridis_c(option = "viridis", na.value = "gray22", begin = 0.6) +
#'    ggplot2::theme_void()
#'  }
#' @export
#'
trees_biomass <- function(
  tree_list
  , crs = NA
  , study_boundary = NA
  , input_landfire_dir = NULL
  , input_foresttype_dir = NULL
  , method = NA
){
  ####################################################################
  # PARSE THE method
  ####################################################################
  # clean method
  method <- dplyr::coalesce(method, "") %>%
    tolower() %>%
    stringr::str_squish()
  # potential methods
  pot_methods <- c("cruz", "landfire") %>% unique()
  find_method <- paste(pot_methods, collapse="|")
  # can i find one?
  which_methods <- stringr::str_extract_all(string = method, pattern = find_method) %>%
    unlist() %>%
    unique()
  # make sure at least one is selected
  # get_list_diff() from get_url_data.R
  n_methods_not <- get_list_diff(pot_methods, which_methods) %>% length()
  if(n_methods_not>=length(pot_methods)){
    stop(paste0(
      "`method` parameter must be one or multiple of:\n"
      , "    "
      , paste(pot_methods, collapse=", ")
    ))
  }

  ####################################################################
  # set up placeholders for return data
  ####################################################################
    stand_cruz <- NULL
    stand_landfire <- NULL

  ####################################################################
  # convert to spatial points data...just do this once
  ####################################################################
    tree_tops <- check_spatial_points(tree_list, crs)

  ####################################################################
  # cruz
  ####################################################################
  if( any(stringr::str_equal(which_methods, "cruz")) ){
    message("attempting to estimate biomass using the `cruz` method")
    # call the function
    safe_trees_biomass_cruz <- purrr::safely(trees_biomass_cruz)
    trees_biomass_cruz_ans <- safe_trees_biomass_cruz(
      tree_list = tree_tops
      , crs = crs
      , study_boundary = study_boundary
      , input_foresttype_dir = input_foresttype_dir
    )
    # check result
    if(is.null(trees_biomass_cruz_ans$error)){
      ############################
      # join new data to tree list
      ############################
      if(inherits(trees_biomass_cruz_ans$result$tree_list, "data.frame")){
        tl_cruz <- trees_biomass_cruz_ans$result$tree_list
        # get names from new df
        names_temp <- c(
          "treeID"
          , get_list_diff(
            names(tl_cruz %>% sf::st_drop_geometry())
            , names(tree_tops %>% sf::st_drop_geometry())
            )
          )
        # join to original data
        tree_tops <- tree_tops %>%
          dplyr::left_join(
            tl_cruz %>%
              sf::st_drop_geometry() %>%
              dplyr::select(dplyr::all_of(names_temp))
            , by = "treeID"
          )
      }
      ############################
      # return raster
      ############################
      if(inherits(trees_biomass_cruz_ans$result$stand_cell_data, "data.frame")){
        stand_cruz <- trees_biomass_cruz_ans$result$stand_cell_data
      }
    }else{
      warning(paste0(
        "Could not get `cruz` biomass estimates:\n"
        , trees_biomass_cruz_ans$error
      ))
    }
  }

  ####################################################################
  # landfire
  ####################################################################
  if( any(stringr::str_equal(which_methods, "landfire")) ){
    message("attempting to estimate biomass using the `landfire` method")
    # call the function
    safe_trees_biomass_landfire <- purrr::safely(trees_biomass_landfire)
    trees_biomass_landfire_ans <- safe_trees_biomass_landfire(
      tree_list = tree_tops
      , crs = crs
      , study_boundary = study_boundary
      , input_landfire_dir = input_landfire_dir
    )
    # check result
    if(is.null(trees_biomass_landfire_ans$error)){
      ############################
      # join new data to tree list
      ############################
      if(inherits(trees_biomass_landfire_ans$result$tree_list, "data.frame")){
        tl_landfire <- trees_biomass_landfire_ans$result$tree_list
        # get names from new df
        names_temp <- c(
          "treeID"
          , get_list_diff(
            names(tl_landfire %>% sf::st_drop_geometry())
            , names(tree_tops %>% sf::st_drop_geometry())
            )
          )
        # join to original data
        tree_tops <- tree_tops %>%
          dplyr::left_join(
            tl_landfire %>%
              sf::st_drop_geometry() %>%
              dplyr::select(dplyr::all_of(names_temp))
            , by = "treeID"
          )
      }
      ############################
      # return raster
      ############################
      if(inherits(trees_biomass_landfire_ans$result$stand_cell_data, "data.frame")){
        stand_landfire <- trees_biomass_landfire_ans$result$stand_cell_data
      }
    }else{
      warning(paste0(
        "Could not get `landfire` biomass estimates:\n"
        , trees_biomass_landfire_ans$error
      ))
    }
  }

  ####################################################################
  # return
  ####################################################################
  return(list(
    tree_list = tree_tops
    , stand_cell_data_landfire = stand_landfire
    , stand_cell_data_cruz = stand_cruz
  ))
}
