#' @title Use outputs from [cloud2trees()] to generate inputs for LANL TREES program
#'
#' @description
#' `cloud2trees_to_lanl_trees()` uses the output from [cloud2trees()] to generate inputs for [LANL TREES](https://github.com/lanl/Trees/)
#' program as a pathway to fire modeling with [Quic-Fire](https://doi.org/10.1016/j.envsoft.2019.104616).
#'
#' The primary input is a directory with outputs from [cloud2trees()]. The default directory written by [cloud2trees()] is `point_cloud_processing_delivery`
#' which must contain (at a minimum):
#'
#' * DTM raster with a name formatted as : "dtm_\*.tif"
#' * Tree list data with the name formatted as : "final_detected_tree_tops_\*.gpkg" (tree points) or "final_detected_crowns_\*.gpkg" (tree crowns)
#' * A study area spatial file that can be read with the `sf` package (see [sf::st_drivers()])
#'
#' @param input_dir directory with outputs from [cloud2trees()]. The default directory written by [cloud2trees()] is `point_cloud_processing_delivery`
#' @param study_boundary sf. The boundary of the study area which is used to determine the outputs
#' @param bbox_aoi logical. Should the study_boundary be transformed to a bounding box instead of it's original shape for determining the objects within the boundary?
#' @param buffer numeric. Buffer to be applied to the study area prior to determining objects within the boundary.
#' Units are determined by the horizontal CRS settings of the tree list data
#' @param topofile character. one of `"flat"` or `"dtm"`:
#'    - "flat" - always flat for QUIC-Fire
#'    - "dtm" - uses the topo.dat file created based on the DTM; potentially for FIRETEC
#' @param cbd_method character. one of `"landfire"` or `"cruz"`:
#' * Tree crown biomass method:
#'    - "landfire" - based on LANDFIRE's Forest Canopy Bulk Density (CBD) data ([trees_biomass_landfire()])
#'    - "cruz" - based on Cruz et al. (2003) canopy fuel stratum equations ([trees_biomass_cruz()])
#' @param output_dir parent directory where new folder `lanl_trees_delivery` will be written for exports
#' @param fuel_litter list. a `list()` or numeric vector `c()`. see default.
#' * must have parameters in order:
#'    - "ilitter" : 0 = no litter, 1 = litter
#'    - "lrho" : litter bulk density
#'    - "lmoisture : litter moisture
#'    - "lss" : litter sizescale
#'    - "ldepth" : litter depth
#' @param fuel_grass list. a `list()` or numeric vector `c()`. see default.
#' * must have parameters in order:
#'    - "igrass" : 0 = no grass, 1 = grass
#'    - "grho" : grass bulk density
#'    - "gmoisture" : grass moisture
#'    - "gss" : grass sizescale
#'    - "gdepth" : grass depth
#'
#' @return Returns a list of objects: tree_list = spatial data frame of individual trees; foresttype_rast = raster of forest types in the area.
#'
#' @export
#'
cloud2trees_to_lanl_trees <- function(
  input_dir
  , study_boundary = NA
  , bbox_aoi = T
  , buffer = 0
  , topofile = "flat"
  , cbd_method = "landfire"
  , output_dir = tempdir()
  # fuel
  , fuel_litter = list(
    "ilitter"  = 0 # 0 = no litter, 1 = litter
    , "lrho"     = 4.667 #litter bulk density
    , "lmoisture"= 0.06 #litter moisture
    , "lss"      = 0.0005 #litter sizescale
    , "ldepth"   = 0.06 #litter depth
  )
  , fuel_grass = list(
    "igrass"   = 0 # 0 = no grass, 1 = grass
    , "grho"     = 1.17 #grass bulk density
    , "gmoisture"= 0.06 #grass moisture
    , "gss"      = 0.0005 #grass sizescale
    , "gdepth"   = 0.27 #grass depth
  )
){
  # potential future parameters
  horizontal_res <- 2
  #######################################
  # checks
  #######################################
    search_dir_final_detected_ans <- search_dir_final_detected(input_dir)
    #### trees
      if(
        is.null(search_dir_final_detected_ans$crowns_flist)
        && is.null(search_dir_final_detected_ans$ttops_flist)
      ){stop(paste0(
        "could not locate tree list data in:\n    "
        , normalizePath(input_dir)
      ))}
    #### dtm
      if(
        is.null(search_dir_final_detected_ans$dtm_flist)
      ){stop(paste0(
        "could not locate DTM raster in:\n    "
        , normalizePath(input_dir)
      ))}
    #### boundary
      if(inherits(study_boundary,"character") && file.exists(study_boundary)){
        study_boundary <- sf::st_read(study_boundary, quiet = T)
      }
      if(
        !inherits(study_boundary,"sf")
        && !inherits(study_boundary,"sfc")
      ){stop("study_boundary must be sf class object")}
      if(is.na(sf::st_crs(study_boundary))){stop("study_boundary does not have a CRS")}
      if(inherits(study_boundary,"sf") && nrow(study_boundary)!=1){
        stop("study_boundary must only have a single record geometry")
      }
      if(inherits(study_boundary,"sfc") && length(study_boundary)!=1){
        stop("study_boundary must only have a single record geometry")
      }
      if(
        !all( sf::st_is(study_boundary, c("POLYGON","MULTIPOLYGON")) )
      ){
        stop("study_boundary must contain POLYGON type geometry only")
      }
    #### topofile
      topofile <- dplyr::coalesce(topofile, "") %>%
        tolower() %>%
        stringr::str_squish()
      # potential methods
      pot_methods <- c("flat", "dtm") %>% unique()
      find_method <- paste(pot_methods, collapse="|")
      # can i find one?
      topofile <- stringr::str_extract_all(string = topofile, pattern = find_method) %>%
        unlist() %>%
        unique()
      # make sure at least one is selected
      # get_list_diff() from get_url_data.R
      n_methods_not <- get_list_diff(pot_methods, topofile) %>% length()
      topofile <- topofile[1] # just get the first
      if(n_methods_not>=length(pot_methods)){
        stop(paste0(
          "`topofile` parameter must be one of:\n"
          , "    "
          , paste(pot_methods, collapse=", ")
        ))
      }
    #### cbd_method
      which_biomass_methods <- check_biomass_method(cbd_method)
      cbd_method <- which_biomass_methods[1]
      if(cbd_method=="landfire"){
        cbd_col_name <- "landfire_tree_kg_per_m3"
      }else if(cbd_method=="cruz"){
        cbd_col_name <- "cruz_tree_kg_per_m3"
      }else{
        stop("could not identify valid `cbd_method`")
      }
    #### fuel
      fuel_litter <- check_fuel_list(fuel_litter, "litter")
      fuel_grass <- check_fuel_list(fuel_grass, "grass")
  #######################################
  # outdir
  #######################################
    outdir <- file.path(normalizePath(output_dir))
    if(!dir.exists(outdir)){
      stop(paste0(
        "could not locate the directory `output_dir` at:\n    "
        , outdir
      ))
    }
    # make the new dir
    outdir <- file.path(outdir,"lanl_trees_delivery")
    if(!dir.exists(outdir)){
      dir.create(outdir, showWarnings = F)
    }else{
      list.files(outdir, recursive = T, full.names = T) %>%
        purrr::map(file.remove)
    }

  #######################################
  # clip_tree_list_aoi
  #######################################
    # get the tree list
    tree_list <- read_trees_flist(input_dir,which_trees = "tops")

    # customize the aoi settings and clip the tree list
    clip_tree_list_aoi_ans <- clip_tree_list_aoi(
      tree_list = tree_list
      , study_boundary = study_boundary
      , bbox_aoi = bbox_aoi
      , buffer = buffer
      , reproject_epsg = NULL
    )
    # clip_tree_list_aoi_ans %>% names()
    # inherits(clip_tree_list_aoi_ans$aoi,"sfc")
    # sf::st_crs(clip_tree_list_aoi_ans$aoi)
    # ggplot() + geom_sf(data=clip_tree_list_aoi_ans$aoi)

  #######################################
  # quicfire_define_domain
  #######################################
    # do the domain thing
    quicfire_domain_df <- quicfire_define_domain(
      sf_data = clip_tree_list_aoi_ans$aoi
      , horizontal_resolution = horizontal_res
      , outdir = outdir
    )
    # quicfire_domain_df
    # quicfire_domain_df$quicfire_domain_df %>% dplyr::glimpse()
    # quicfire_domain_df$quicfire_domain_df %>% sf::st_crs(parameters = T)
    # nrow(quicfire_domain_df$quicfire_domain_df)

  #######################################
  # quicfire_dtm_topofile
  #######################################
    # don't forget the fortran flipped over rast
    quicfire_dtm_topofile_ans <- quicfire_dtm_topofile(
      dtm_rast = search_dir_final_detected_ans$dtm_flist[1]
      , horizontal_resolution = horizontal_res
      , study_boundary = quicfire_domain_df$quicfire_domain_df
      , outdir = outdir
    )
    # quicfire_dtm_topofile_ans
    # quicfire_dtm_topofile_ans$topofile_path
    # terra::plot(quicfire_dtm_topofile_ans$dtm)
    # terra::plot(
    #   quicfire_domain_df$quicfire_domain_df %>%
    #     sf::st_transform(terra::crs(quicfire_dtm_topofile_ans$dtm)) %>%
    #     terra::vect()
    #   , add = T, border = "blue", col = NA
    #   , lwd = 22
    # )

  #######################################
  # make_lanl_trees_input
  #######################################
    if(topofile=="dtm"){
      topofile <- quicfire_dtm_topofile_ans$topofile_path
    }else{
      topofile <- "flat"
    }

    make_lanl_trees_input_ans <- make_lanl_trees_input(
      tree_list = clip_tree_list_aoi_ans$tree_list
      , quicfire_domain_df = quicfire_domain_df$quicfire_domain_df
      , topofile = topofile
      , cbd_col_name = cbd_col_name # "cruz_tree_kg_per_m3" "landfire_tree_kg_per_m3"
      , horizontal_resolution = horizontal_res
      , outdir = outdir
      # fuel settings
      , fuel_litter = fuel_litter
      , fuel_grass = fuel_grass
    )

}
