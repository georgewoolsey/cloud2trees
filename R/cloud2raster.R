#' @title Use raw .las|.laz files to generate CHM, DTM, and Normalized .las files
#'
#' @description
#' `cloud2raster()` is an all-in-one function to process raw .las|.laz files
#' to generate a CHM raster (.tif), a DTM raster (.tif), and .las files which have been height normalized.
#' The order of operations is:
#'
#' * Tile the raw point cloud to work with smaller chunks and reduce the potential for memory issues with high density clouds using [chunk_las_catalog()]
#' * Classify the point cloud using [lasR::classify_with_csf()]
#' * Remove outlier points using [lasR::classify_with_ivf()]
#' * Produce a triangulation of the ground points (meshed DTM) using [lasR::triangulate()]
#' * Rasterize the result of the Delaunay triangulation using [lasR::rasterize()]
#' * Height normalize the point cloud using either the DTM or the triangulation [lasR::transform_with()]
#' * Use the height normalized point cloud to create the CHM based on the highest point in a pixel using [lasR::rasterize()]
#' * Pits and spikes filling of the CHM raster using [lasR::pit_fill()]
#' * Smooth the CHM raster tile gaps using [terra::focal()]
#'
#' @inheritParams create_project_structure
#' @inheritParams chunk_las_catalog
#' @param keep_intrmdt logical. this process writes intermediate data to the disk, keep those intermediate files (classfied, normalized, stem las files)?
#' @param dtm_res_m numeric. The desired resolution of the DTM produced in meters.
#' @param chm_res_m numeric. The desired resolution of the CHM produced in meters.
#' @param min_height numeric. Set the minimum height (m) for individual tree detection
#' @param max_height numeric. Set the maximum height (m) for the canopy height model
#'
#' @references
#' https://r-lidar.github.io/lasR/index.html
#' https://r-lidar.github.io/lidRbook/norm.html
#'
#'
#' @return Returns the goods.
#' Exports files of the goods to new folders "point_cloud_processing_delivery" and "point_cloud_processing_temp" in the
#' `output_dir` defined by the user in the function call.
#'
#' @examples
#'  \dontrun{
#'  f <- "../lasdata"
#'  chunk_las_catalog(folder = f, outfolder = getwd())
#'  }
#' @export
#'
cloud2raster <- function(
  output_dir
  , input_las_dir
  , input_treemap_dir = paste0(system.file(package = "cloud2trees"),"/extdata/treemap")
  , accuracy_level = 2
  , max_ctg_pts = 70e6
  , max_area_m2 = 90e6
  , transform = FALSE
  , new_crs = NA
  , old_crs = NA
  , keep_intrmdt = F
  , dtm_res_m = 1
  , chm_res_m = 0.25
  , min_height = 2
  , max_height = 70
){
  ######################################
  # execute pipeline
  ######################################
    # set lasR parallel processing options as of lasR 0.4.2
      lasR::set_parallel_strategy(
        strategy = lasR::concurrent_points(
            ncores = max(lasR::ncores()-1, lasR::half_cores())
          )
      )

    # create safe function to capture error and map over
      safe_lasr_pipeline = purrr::safely(lasr_pipeline)

    # map over processing grids
      lasr_ans_list =
        process_data$processing_grid %>%
        unique() %>%
        purrr::map(safe_lasr_pipeline)

    # check for errors other than too few points for triangulation which happens on edge chunks with few points and also those with copious noise
      error_list_temp =
        lasr_ans_list %>%
        purrr::transpose() %>%
        purrr::pluck("error") %>%
        unlist() %>%
        purrr::pluck("message")

      if(length(error_list_temp)>0){
        has_errors_temp = error_list_temp %>%
          dplyr::tibble() %>%
          dplyr::rename(message=1) %>%
          dplyr::mutate(
            is_tri_error = dplyr::case_when(
              stringr::str_detect(tolower(message), "impossible") &
                stringr::str_detect(tolower(message), "triangulation") ~ 1
              , T ~ 0
            )
            , sum_is_tri_error = sum(is_tri_error)
            , pct_tri = sum_is_tri_error/length(unique(process_data$processing_grid))
            , keep_it = dplyr::case_when(
              is_tri_error==1 & pct_tri<0.5 ~ F
              , T ~ T
            )
          ) %>%
          dplyr::filter(keep_it==T)
      }else{has_errors_temp = dplyr::tibble()}

    if(nrow(has_errors_temp)>0){stop("lasR processing failed due to:\n", has_errors_temp$message[1])}

    # clean up
      remove(list = ls()[grep("_temp",ls())])
      gc()
  ###################################
  # lasR cleanup and polish
  ###################################
    ###############
    # DTM
    ###############
      # output name
      dtm_file_name = paste0(config$delivery_dir, "/dtm_", desired_dtm_res, "m.tif")
      # read
      rast_list_temp = list.files(config$dtm_dir, pattern = ".*\\.(tif|tiff)$", full.names = T) %>% purrr::map(function(x){terra::rast(x)})
      # mosaic
      dtm_rast = terra::sprc(rast_list_temp) %>% terra::mosaic(fun = "mean")
      # # fill empty cells
      #   #### this is not needed anymore as empty cells resulting from rasterizing with a triangulation are fixed
      #   #### see: https://github.com/r-lidar/lasR/issues/18
      #   dtm_rast = dtm_rast %>%
      #     terra::crop(
      #       las_ctg@data$geometry %>%
      #         sf::st_union() %>%
      #         terra::vect() %>%
      #         terra::project(terra::crs(dtm_rast))
      #     ) %>%
      #     terra::mask(
      #       las_ctg@data$geometry %>%
      #         sf::st_union() %>%
      #         terra::vect() %>%
      #         terra::project(terra::crs(dtm_rast))
      #     ) %>%
      #     terra::focal(
      #       w = 3
      #       , fun = "mean"
      #       , na.rm = T
      #       # na.policy Must be one of:
      #         # "all" (compute for all cells)
      #         # , "only" (only for cells that are NA)
      #         # , or "omit" (skip cells that are NA).
      #       , na.policy = "only"
      #     )
      #   # dtm_rast %>% terra::crs()
      #   # dtm_rast %>% terra::plot()
      # set crs
        terra::crs(dtm_rast) = proj_crs
      # write to delivery directory
        terra::writeRaster(
          dtm_rast
          , filename = dtm_file_name
          , overwrite = T
        )
    ###############
    # chm
    ###############
      # output name
      chm_file_name = paste0(config$delivery_dir, "/chm_", desired_chm_res, "m.tif")
      # read
      rast_list_temp = list.files(config$chm_dir, pattern = ".*\\.(tif|tiff)$", full.names = T) %>% purrr::map(function(x){terra::rast(x)})
      # mosaic
      chm_rast = terra::sprc(rast_list_temp) %>% terra::mosaic(fun = "max")
      # set crs
      terra::crs(chm_rast) = proj_crs
      # fill empty cells
        # this helps to smooth out tile gaps which leads to too many trees being detected during itd
        chm_rast = chm_rast %>%
          terra::crop(
            las_ctg@data$geometry %>%
              sf::st_union() %>%
              terra::vect()
          ) %>%
          terra::mask(
            las_ctg@data$geometry %>%
              sf::st_union() %>%
              terra::vect()
          ) %>%
          terra::focal(
            w = 3
            , fun = "mean"
            , na.rm = T
            # na.policy Must be one of:
              # "all" (compute for all cells)
              # , "only" (only for cells that are NA)
              # , or "omit" (skip cells that are NA).
            , na.policy = "only"
          )
        # chm_rast %>% terra::crs()
        # chm_rast %>% terra::plot()

      # write to delivery directory
        terra::writeRaster(
          chm_rast
          , filename = chm_file_name
          , overwrite = T
        )

    # create spatial index files (.lax)
      # classify
      create_lax_for_tiles(
        las_file_list = list.files(config$las_classify_dir, pattern = ".*\\.(laz|las)$", full.names = T)
      )
      # normalize
      normalize_flist = create_lax_for_tiles(
        las_file_list = list.files(config$las_normalize_dir, pattern = ".*\\.(laz|las)$", full.names = T)
      )

    # clean up
      remove(list = ls()[grep("_temp",ls())])
      gc()

    # # # plots
    # dtm_rast %>%
    #   # terra::aggregate(fact = 4) %>%
    #   as.data.frame(xy = T) %>%
    #   dplyr::rename(f=3) %>%
    #   ggplot() +
    #     geom_tile(aes(x=x,y=y,fill=f)) +
    #     geom_sf(data = las_ctg@data$geometry, fill = NA) +
    #     scale_fill_viridis_c(option = "viridis") +
    #     labs(fill = "meters") +
    #     theme_void()
    # chm_rast %>%
    #   as.data.frame(xy = T) %>%
    #   dplyr::rename(f=3) %>%
    #   ggplot() +
    #     geom_tile(aes(x=x,y=y,fill=f)) +
    #     geom_sf(data = las_ctg@data$geometry, fill = NA) +
    #     scale_fill_viridis_c(option = "plasma") +
    #     labs(fill = "meters") +
    #     theme_void()
}
