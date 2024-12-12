#' Create project structure
#'
#' @description
#' Function to generate nested project directories based on the user-defined directory to create output file structure
#'
#' @param output_dir parent directory where new folders `point_cloud_processing_delivery` and `point_cloud_processing_temp` will be written for exports
#' @param input_las_dir directory where .las|.laz point cloud data exists...program will search all sub-directories for all .las|.laz files and process them as one
#' @inheritParams find_ext_data
#'
#' @return A data.frame.
#'
#' @keywords internal
#'
create_project_structure <- function(
    output_dir
    , input_las_dir
    , input_treemap_dir = NULL
    , input_foresttype_dir = NULL
){
  # check if folder contains las files directly
    chk <- input_las_dir %>%
      tolower() %>%
      basename() %>%
      stringr::str_detect(".*\\.(laz|las)$")
    # if it's a file or file list...get the directory
    if(max(chk)==1){
      input <- input_las_dir %>%
      stringr::str_subset(".*\\.(laz|las)$") %>%
      dirname() %>%
      normalizePath()
    }else{ # otherwise just list the folder
      input <- normalizePath(input_las_dir)
    }

  ###___________________________________________________###
  ### Set output delivery directory
  ###___________________________________________________###
  temp_dir <- file.path(normalizePath(output_dir), "point_cloud_processing_temp")
  delivery_dir <- file.path(normalizePath(output_dir), "point_cloud_processing_delivery")
  ### set output directory for temporary files
  las_grid_dir <- file.path(temp_dir, "00_grid")
  las_classify_dir <- file.path(temp_dir, "01_classify")
  las_normalize_dir <- file.path(temp_dir, "02_normalize")
  dtm_dir <- file.path(temp_dir, "03_dtm")
  chm_dir <- file.path(temp_dir, "04_chm")
  treels_dir <- file.path(temp_dir, "05_treels")
  ### Create the directories
  dir.create(delivery_dir, showWarnings = FALSE)
  dir.create(temp_dir, showWarnings = FALSE)
  dir.create(las_grid_dir, showWarnings = FALSE)
  dir.create(las_classify_dir, showWarnings = FALSE)
  dir.create(las_normalize_dir, showWarnings = FALSE)
  dir.create(dtm_dir, showWarnings = FALSE)
  dir.create(chm_dir, showWarnings = FALSE)
  dir.create(treels_dir, showWarnings = FALSE)
  ###______________________________###
  ### create return ###
  ###______________________________###
  config <- dplyr::tibble(
    input_las_dir = input
    , input_treemap_dir = dplyr::coalesce(input_treemap_dir, as.character(NA))
    , input_foresttype_dir = dplyr::coalesce(input_foresttype_dir, as.character(NA))
    , delivery_dir = delivery_dir
    , temp_dir = temp_dir
    , las_grid_dir = las_grid_dir
    , las_classify_dir = las_classify_dir
    , las_normalize_dir = las_normalize_dir
    , dtm_dir = dtm_dir
    , chm_dir = chm_dir
    , treels_dir = treels_dir
    , is_input_file_list = max(chk)==1
  )

  #config
  # remove all files in delivery and temp
  list.files(config$temp_dir, recursive = T, full.names = T) %>%
      purrr::map(file.remove)
  list.files(config$delivery_dir, recursive = T, full.names = T) %>%
      purrr::map(file.remove)

  ### Return config
  return(config)

}
