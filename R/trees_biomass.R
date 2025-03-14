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
#' # use the tree list that ships with the package
#' f <- system.file(package = "cloud2trees", "extdata", "crowns_poly.gpkg")
#' tl <- sf::st_read(f)
#' tl %>% dplyr::glimpse()
#' # call trees_biomass and get multiple biomass estimates
#' trees_biomass_ans <- trees_biomass(tree_list = tl, method = c("landfire","cruz"))
#' # what did we get back?
#' trees_biomass_ans %>% names()
#' # check out the tree list
#' trees_biomass_ans$tree_list %>% dplyr::glimpse()
#' # check out the landfire stand data
#' trees_biomass_ans$stand_cell_data_landfire %>% dplyr::filter(trees>0) %>% dplyr::glimpse()
#' # plot tree landfire crown biomass estimate
#' trees_biomass_ans$tree_list %>%
#'   ggplot2::ggplot(
#'     mapping = ggplot2::aes(
#'       x = tree_height_m
#'       , y = landfire_crown_biomass_kg
#'       , color = crown_area_m2
#'     )
#'   ) +
#'   ggplot2::geom_point()
#' # plot tree cruz crown biomass estimate
#' trees_biomass_ans$tree_list %>%
#'   ggplot2::ggplot(
#'     mapping = ggplot2::aes(
#'       x = tree_height_m
#'       , y = cruz_crown_biomass_kg
#'       , color = crown_area_m2
#'     )
#'   ) +
#'   ggplot2::geom_point()
#' # plot tree landfire vs. cruz crown biomass estimate
#' trees_biomass_ans$tree_list %>%
#'   ggplot2::ggplot(
#'     mapping = ggplot2::aes(
#'       x = landfire_crown_biomass_kg, y = cruz_crown_biomass_kg
#'     )
#'   ) +
#'   ggplot2::geom_abline(lwd = 1.5) +
#'   ggplot2::geom_smooth(method = "lm", se=F, color = "gray", linetype = "dashed") +
#'   ggplot2::geom_point(ggplot2::aes(color = tree_height_m)) +
#'   ggplot2::scale_x_continuous(
#'     limits = c(0
#'       , max(trees_biomass_ans$tree_list$cruz_crown_biomass_kg)
#'     )
#'   ) +
#'   ggplot2::scale_y_continuous(
#'     limits = c(0
#'       , max(trees_biomass_ans$tree_list$cruz_crown_biomass_kg)
#'     )
#'   )
#' # get the projection for the stand cell data
#' epsg_code <- trees_biomass_ans$stand_cell_data_landfire$rast_epsg_code[1] %>% as.numeric()
#' # plot the stand cell data with trees overlaid
#' trees_biomass_ans$stand_cell_data_landfire %>%
#'   dplyr::filter(trees>0) %>%
#'   ggplot2::ggplot() +
#'   ggplot2::geom_tile(ggplot2::aes(x=x,y=y,fill = landfire_stand_kg_per_m3), color = "gray44") +
#'   ggplot2::geom_text(ggplot2::aes(x=x,y=y,label = trees), color = "white") +
#'   ggplot2::geom_sf(
#'     data = trees_biomass_ans$tree_list %>% sf::st_transform(crs = epsg_code)
#'     , ggplot2::aes(color = cruz_crown_biomass_kg)
#'   ) +
#'   ggplot2::labs(
#'     fill="stand kg/m3", color = "landfire\ncrown kg"
#'     , caption = "# trees shown in cell"
#'   ) +
#'   ggplot2::scale_fill_viridis_c(option = "rocket", na.value = "gray", direction = -1) +
#'   ggplot2::scale_color_viridis_c(option = "viridis", na.value = "gray22", begin = 0.6) +
#'   ggplot2::theme_void()
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
  , max_crown_kg_per_m3 = 2
){
  ####################################################################
  # PARSE THE method
  # will throw error if not in available options
  ####################################################################
    which_biomass_methods <- check_biomass_method(method)

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
  if( any(stringr::str_equal(which_biomass_methods, "cruz")) ){
    message("attempting to estimate biomass using the `cruz` method")
    # call the function
    safe_trees_biomass_cruz <- purrr::safely(trees_biomass_cruz)
    trees_biomass_cruz_ans <- safe_trees_biomass_cruz(
      tree_list = tree_tops
      , crs = crs
      , study_boundary = study_boundary
      , input_foresttype_dir = input_foresttype_dir
      , max_crown_kg_per_m3 = max_crown_kg_per_m3
    )
    # check result
    if(is.null(trees_biomass_cruz_ans$error)){
      ############################
      # join new data to tree list
      ############################
      if(inherits(trees_biomass_cruz_ans$result$tree_list, "data.frame")){
        tl_cruz <- trees_biomass_cruz_ans$result$tree_list
        # remove columns created in trees_biomass_cruz if already existed in orig data
        tree_tops <- tree_tops %>%
          clean_biomass_cols(method = "cruz")
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
  if( any(stringr::str_equal(which_biomass_methods, "landfire")) ){
    message("attempting to estimate biomass using the `landfire` method")
    # call the function
    safe_trees_biomass_landfire <- purrr::safely(trees_biomass_landfire)
    trees_biomass_landfire_ans <- safe_trees_biomass_landfire(
      tree_list = tree_tops
      , crs = crs
      , study_boundary = study_boundary
      , input_landfire_dir = input_landfire_dir
      , max_crown_kg_per_m3 = max_crown_kg_per_m3
    )
    # check result
    if(is.null(trees_biomass_landfire_ans$error)){
      ############################
      # join new data to tree list
      ############################
      if(inherits(trees_biomass_landfire_ans$result$tree_list, "data.frame")){
        tl_landfire <- trees_biomass_landfire_ans$result$tree_list
        # remove columns created in trees_biomass_landfire if already existed in orig data
        tree_tops <- tree_tops %>%
          clean_biomass_cols(method = "landfire")
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
