#' create project structure
#' @description
#' Function to generate nested project directories based on the user-defined directory to create output file structure
#'
#' @param output_dir parent directory where new folders `point_cloud_processing_delivery` and `point_cloud_processing_temp` will be written for exports
#' @param input_las_dir directory where .las|.laz point cloud data exists...program will search all sub-directories for all .las|.laz files and process them as one
#' @param input_treemap_dir directory where Treemap 2016 exists. Use `get_treemap()` first.
#'
#' @return A data.frame.
#'
#' @keywords internal
#'
create_project_structure = function(
    output_dir, input_las_dir,
    input_treemap_dir = paste0(system.file(package = "cloud2trees"),"/extdata/treemap")
){
  # check treemap data... see get_treemap()
  f <- toupper(list.files(input_treemap_dir))
  if(length(f)==0){f <- ""}
  if(
    max(grepl("TREEMAP2016.TIF", f))==0 | max(grepl("TREEMAP2016_TREE_TABLE.CSV", f))==0
  ){
    stop("Treemap data has not been downloaded to package contents. Use `get_treemap()` first.")
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
  las_stem_norm_dir <- file.path(temp_dir, "05_stem_norm")
  las_stem_dir <- file.path(temp_dir, "06_las_stem")
  stem_poly_tile_dir <- file.path(temp_dir, "07_stem_poly_tile")
  ### Create the directories
  dir.create(delivery_dir, showWarnings = FALSE)
  dir.create(temp_dir, showWarnings = FALSE)
  dir.create(las_grid_dir, showWarnings = FALSE)
  dir.create(las_classify_dir, showWarnings = FALSE)
  dir.create(las_normalize_dir, showWarnings = FALSE)
  dir.create(dtm_dir, showWarnings = FALSE)
  dir.create(chm_dir, showWarnings = FALSE)
  dir.create(las_stem_norm_dir, showWarnings = FALSE)
  dir.create(las_stem_dir, showWarnings = FALSE)
  dir.create(stem_poly_tile_dir, showWarnings = FALSE)
  ###______________________________###
  ### Set names of the directories ###
  ###______________________________###
  names(output_dir) <- "output_dir"
  names(input_las_dir) <- "input_las_dir"
  names(input_treemap_dir) <- "input_treemap_dir"
  names(delivery_dir) <- "delivery_dir"
  names(temp_dir) <- "temp_dir"
  names(las_grid_dir) <- "las_grid_dir"
  names(las_classify_dir) <- "las_classify_dir"
  names(las_normalize_dir) <- "las_normalize_dir"
  names(dtm_dir) <- "dtm_dir"
  names(chm_dir) <- "chm_dir"
  names(las_stem_norm_dir) <- "las_stem_norm_dir"
  names(las_stem_dir) <- "las_stem_dir"
  names(stem_poly_tile_dir) <- "stem_poly_tile_dir"

  ###______________________________###
  ### Append to output config list ###
  ###______________________________###

  config <- cbind(
    output_dir, input_las_dir, input_treemap_dir
    , delivery_dir, temp_dir, las_grid_dir, las_classify_dir
    , las_normalize_dir, dtm_dir, chm_dir
    , las_stem_norm_dir, las_stem_dir, stem_poly_tile_dir
  )

  config <- as.data.frame(config)

  #config
  # remove all files in delivery and temp
  list.files(config$temp_dir, recursive = T, full.names = T) %>%
      purrr::map(file.remove)
  list.files(config$delivery_dir, recursive = T, full.names = T) %>%
      purrr::map(file.remove)

  ### Return config
  return(config)

}
