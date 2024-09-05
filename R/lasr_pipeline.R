#' Create lasR pipeline to process a las grid tile created via `chunk_las_catalog`
#' @description
#' Create lasR pipeline to process a las grid tile created via `chunk_las_catalog`
#'
#' @param processing_grid_num numeric. processing_grid column in the data.frame created via `chunk_las_catalog`
#' @param process_data dataframe. data.frame created via `chunk_las_catalog`
#' @param keep_intrmdt logical. this process writes intermediate data to the disk, keep those intermediate files (classfied, normalized, stem las files)?
#' @param dtm_res_m numeric. The desired resolution of the DTM produced in meters.
#' @param chm_res_m numeric. The desired resolution of the CHM produced in meters.
#' @param min_height numeric. Set the minimum height (m) for individual tree detection
#' @param max_height numeric. Set the maximum height (m) for the canopy height model
#' @param dtm_dir string. The path of a folder to write the tiled DTM files to.
#' @param chm_dir string. The path of a folder to write the tiled CHM files to.
#' @param classify_dir string. The path of a folder to write the classified .las files to.
#' @param normalize_dir string. The path of a folder to write the normalized .las files to.
#'
#' @references
#' https://r-lidar.github.io/lasR/index.html
#'
#' @return A `lasR` pipeline
#'
#' @keywords internal
#'
lasr_pipeline <- function(
  processing_grid_num = 1 # will map over this
  # these parameters are meant to be set once
  , process_data
  , keep_intrmdt = F
  , dtm_res_m = 1
  , chm_res_m = 0.25
  , min_height = 2
  , max_height = 70
  , dtm_dir = getwd()
  , chm_dir = getwd()
  , classify_dir = getwd()
  , normalize_dir = getwd()
){
  #################################################
  # setup to pass to lasR functions
  #################################################
    # output files
    dtm_file_name <- paste0(normalizePath(dtm_dir), "/", processing_grid_num,"_dtm_", dtm_res_m, "m.tif")
    chm_file_name <- paste0(normalizePath(chm_dir), "/", processing_grid_num,"_chm_", chm_res_m, "m.tif")
    # files to process
    flist <- process_data %>%
      dplyr::filter(processing_grid == processing_grid_num) %>%
      dplyr::pull(filename)
    # set up catalog
    tile_las_ctg <- lidR::readLAScatalog(flist)
    lidR::opt_chunk_buffer(tile_las_ctg) <- 10
    lidR::opt_chunk_size(tile_las_ctg) <- 0
    # fraction for triangulation
    frac_for_tri <- process_data %>%
      dplyr::filter(processing_grid == processing_grid_num) %>%
      dplyr::pull(pts_m2_factor_ctg) %>%
      .[1]
    # pull accuracy info from data.frame
    normalization_accuracy <- process_data %>%
      dplyr::filter(processing_grid == processing_grid_num) %>%
      dplyr::pull(normalization_accuracy) %>%
      .[1]
    # !!! testing with low point density lidar data revealed that the default pit_fill algorithm was too aggressive
    # ... decrease the lasR::pit_fill Size of the Laplacian filter kernel (integer value, in pixels) for low density point clouds
    lap_size_pts <- process_data %>%
      dplyr::filter(processing_grid == processing_grid_num) %>%
      dplyr::pull(pts_m2) %>%
      .[1]
  #################################################
  # DEFINE ALL lasR pipeline steps
  #################################################
    ###################
    # read with filter
    ###################
      lasr_read <- lasR::reader_las(filter = "-drop_noise -drop_duplicates")
    ###################
    # classify
    ###################
      # Classify ground points
      # library(RCSF) # for the Cloth simulation filtering (CSF)
      # (Zhang et al 2016) algorithm to classify points
      lasr_classify <- lasR::classify_with_csf(
        slope_smooth = FALSE
        , class_threshold = 0.5
        , cloth_resolution = 0.5
        , rigidness = 1L
        , iterations = 500L
        , time_step = 0.65
        , class = 2L
        , filter = "-drop_noise -drop_duplicates"
      )
    ###################
    # denoise
    ###################
      # classify isolated points
      lasr_denoise <- lasR::classify_with_ivf(res =  5, n = 6)
    ###################
    # write
    ###################
      # write class
      lasr_write_classify <- lasR::write_las(
        ofile = paste0(normalizePath(classify_dir), "/*_classify.las")
        , filter = "-drop_noise"
        , keep_buffer = F
      )
      # write norm
      lasr_write_normalize <- lasR::write_las(
        filter = "-drop_z_below 0 -drop_class 18"
        , ofile = paste0(normalizePath(normalize_dir), "/*_normalize.las")
        , keep_buffer = F
      )
  #################################################
  # build lasR pipeline
  #################################################
  if(keep_intrmdt == T){
    # pipeline
      lasr_pipeline_temp <- lasr_read +
        lasr_classify +
        lasr_denoise +
        lasr_write_classify +
        lasr_dtm_norm(
          dtm_file_name = dtm_file_name
          , frac_for_tri = frac_for_tri
          , dtm_res = dtm_res_m
          , norm_accuracy = normalization_accuracy
        ) +
        lasr_write_normalize +
        lasr_chm(
          chm_file_name = chm_file_name
          , chm_res = chm_res_m
          , min_height_m = min_height
          , max_height_m = max_height
          , lap_sz = ifelse(lap_size_pts<20,2,3)
        )
  }else{
    # pipeline
      lasr_pipeline_temp <- lasr_read +
        lasr_classify +
        lasr_denoise +
        lasr_dtm_norm(
          dtm_file_name = dtm_file_name
          , frac_for_tri = frac_for_tri
          , dtm_res = dtm_res_m
          , norm_accuracy = normalization_accuracy
        ) +
        lasr_write_normalize +
        lasr_chm(
          chm_file_name = chm_file_name
          , chm_res = chm_res_m
          , min_height_m = min_height
          , max_height_m = max_height
          , lap_sz = ifelse(lap_size_pts<20,2,3)
        )
  }
  # message
  message(
    "starting lasR processing of processing grid "
    , processing_grid_num
    , " ("
    , process_data %>%
      dplyr::filter(processing_grid == processing_grid_num) %>%
      dplyr::pull(processing_grid_tot_pts) %>%
      .[1] %>%
      scales::comma(accuracy = 1)
    , " pts) at ..."
    , Sys.time()
  )
  # lasR execute
  lasr_ans <- lasR::exec(
    lasr_pipeline_temp
    , on = tile_las_ctg
    , with = list(
      progress = T
    )
  )
  # sleep
  Sys.sleep(3) # sleep to give c++ stuff time to reset
  return(lasr_ans)
}
