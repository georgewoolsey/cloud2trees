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
#' * Rasterize the result of the Delaunay triangulation using [lasR::rasterize()] to create a DTM
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
#' @param overwrite logical. Should the output files in the `point_cloud_processing_delivery` directory from previous iterations be deleted?
#'
#' @references
#' [https://r-lidar.github.io/lasR/index.html](https://r-lidar.github.io/lasR/index.html)
#' [https://r-lidar.github.io/lidRbook/normalization.html](https://r-lidar.github.io/lidRbook/normalization.html)
#'
#'
#' @return Returns the goods.
#' Exports files of the goods to new folders "point_cloud_processing_delivery" and "point_cloud_processing_temp" in the
#' `output_dir` defined by the user in the function call.
#'
#' @examples
#'  \dontrun{
#'  # test las file but this could also be a directory path with >1 .las|.laz files
#'  i <- system.file("extdata", "MixedConifer.laz", package="lidR")
#'  # run it
#'  r <- cloud2trees::cloud2raster(output_dir = tempdir(), input_las_dir = i)
#'  # what is it?
#'  r %>% names()
#'  # there's a DTM
#'  r$dtm_rast %>% terra::plot()
#'  # there's a CHM
#'  r$chm_rast %>% terra::plot()
#'  # there's a data.frame with the file structure for the project
#'  r$create_project_structure_ans %>% dplyr::glimpse()
#'  # there's a information detailing how the point cloud was processed
#'  r$chunk_las_catalog_ans$process_data %>% dplyr::glimpse()
#'  r$chunk_las_catalog_ans$is_chunked_grid
#'  r$chunk_las_catalog_ans$las_ctg@data %>% dplyr::glimpse()
#'  # there's a list of the height normalized .las files created
#'  r$normalize_flist
#'  }
#' @export
#'
cloud2raster <- function(
  output_dir
  , input_las_dir
  , input_treemap_dir = NULL
  , input_foresttype_dir = NULL
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
  , overwrite = TRUE
){
  ######################################
  # Configure File Structure
  ######################################
  # find external data
  find_ext_data_ans <- find_ext_data(
    input_treemap_dir = input_treemap_dir
    , input_foresttype_dir = input_foresttype_dir
  )

  config <- create_project_structure(
    output_dir = output_dir
    , input_las_dir = input_las_dir
    , input_treemap_dir = find_ext_data_ans$treemap_dir
    , input_foresttype_dir = find_ext_data_ans$foresttype_dir
  )

  # remove all files in delivery and temp
  list.files(config$temp_dir, recursive = T, full.names = T) %>%
      purrr::map(file.remove)
  if(overwrite == TRUE){
    list.files(config$delivery_dir, recursive = T, full.names = T) %>%
        purrr::map(file.remove)
  }

  ######################################
  # Tile raw las files to work with smaller chunks
  ######################################
  chunk_las_catalog_ans <- chunk_las_catalog(
    folder = ifelse(
      config$is_input_file_list == T
      , input_las_dir
      , config$input_las_dir
    )
    , outfolder = config$las_grid_dir
    , accuracy_level = accuracy_level
    , max_ctg_pts = max_ctg_pts
    , max_area_m2 = max_area_m2
    , transform = transform
    , new_crs = new_crs
    , old_crs = old_crs
  )

  ######################################
  # execute lasR pipeline
  ######################################
    # set lasR parallel processing options as of lasR 0.10.1
      lasR::set_parallel_strategy(
        strategy = lasR::concurrent_points(
            ncores = max(lasR::ncores()-2, lasR::half_cores())
          )
      )

    # create safe function to capture error and map over
      # lasr_pipeline <- lasr_pipeline() # is this needed?
      safe_lasr_pipeline <- purrr::safely(lasr_pipeline)

    # map over processing grids
      lasr_ans_list <-
        chunk_las_catalog_ans$process_data$processing_grid %>%
        unique() %>%
        purrr::map(\(x) safe_lasr_pipeline(
          processing_grid_num = x
          , process_data = chunk_las_catalog_ans$process_data
          , keep_intrmdt = keep_intrmdt
          , dtm_res_m = dtm_res_m
          , chm_res_m = chm_res_m
          , min_height = min_height
          , max_height = max_height
          , dtm_dir = config$dtm_dir
          , chm_dir = config$chm_dir
          , classify_dir = config$las_classify_dir
          , normalize_dir = config$las_normalize_dir
        ))

    # check for errors other than too few points for triangulation which happens on edge chunks with few points and also those with copious noise
      error_list_temp <-
        lasr_ans_list %>%
        purrr::transpose() %>%
        purrr::pluck("error") %>%
        unlist() %>%
        purrr::pluck("message")

      if(length(error_list_temp)>0){
        has_errors_temp <- error_list_temp %>%
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
      }else{has_errors_temp <- dplyr::tibble()}

    if(nrow(has_errors_temp)>0){stop("lasR processing failed due to:\n", has_errors_temp$message[1])}

  ###################################
  # lasR cleanup and polish
  ###################################
    ###############
    # DTM
    ###############
      # output name
      dtm_file_name <- paste0(config$delivery_dir, "/dtm_", dtm_res_m, "m.tif")
      # read
      rast_list_temp <- list.files(config$dtm_dir, pattern = ".*\\.(tif|tiff)$", full.names = T) %>%
        purrr::map(function(x){terra::rast(x)})
      # mosaic
      dtm_rast <- terra::sprc(rast_list_temp) %>% terra::mosaic(fun = "mean")
      # set crs
        terra::crs(dtm_rast) <- chunk_las_catalog_ans$process_data$proj_crs[1]
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
      chm_file_name <- paste0(config$delivery_dir, "/chm_", chm_res_m, "m.tif")
      # read
      rast_list_temp <- list.files(config$chm_dir, pattern = ".*\\.(tif|tiff)$", full.names = T) %>%
        purrr::map(function(x){terra::rast(x)})
      # mosaic
      chm_rast <- terra::sprc(rast_list_temp) %>% terra::mosaic(fun = "max")
      # set crs
      terra::crs(chm_rast) <- chunk_las_catalog_ans$process_data$proj_crs[1]
      # fill empty cells
        # this helps to smooth out tile gaps which leads to too many trees being detected during itd
        chm_rast <- chm_rast %>%
          terra::crop(
            chunk_las_catalog_ans$las_ctg@data$geometry %>%
              sf::st_union() %>%
              terra::vect()
          ) %>%
          terra::mask(
            chunk_las_catalog_ans$las_ctg@data$geometry %>%
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
      normalize_flist <- create_lax_for_tiles(
        las_file_list = list.files(config$las_normalize_dir, pattern = ".*\\.(laz|las)$", full.names = T)
      )

    # return
    return(list(
      dtm_rast = terra::rast(dtm_file_name)
      , chm_rast = terra::rast(chm_file_name)
      , create_project_structure_ans = config
      , chunk_las_catalog_ans = chunk_las_catalog_ans
      , normalize_flist = normalize_flist
    ))
}
