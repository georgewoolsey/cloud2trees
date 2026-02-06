#' @title Workflow for slash pile detection, quantification, and spectral validation
#'
#' @description A master wrapper function that executes the full slash pile detection
#' and quantification framework. It identifies structural candidates using [piles_detect()] and
#' subsequently applies the spectral threshold voting system via
#' [piles_spectral_filter()] to generate final pile predictions.
#'
#' @inheritParams piles_detect
#' @param rgb_rast A multi-band SpatRaster containing RGB imagery.
#' @param red_band_idx Integer. Index of the red band.
#' @param green_band_idx Integer. Index of the green band.
#' @param blue_band_idx Integer. Index of the blue band.
#' @param spectral_weight Integer (0 to 6). The consensus threshold
#' requiring a specific number of index thresholds to be met for a candidate
#' pile to be retained.
#' @param filter_return Logical. If TRUE (default), the function returns
#' only polygons where the "inrange_th_votes" is greater than or equal to
#' the `spectral_weight`. If FALSE, the function returns the entire input
#' dataset with columns for each calculated spectral index and the
#' "inrange_th_votes" column for post-processing and analysis.
#'
#' @return If `outfile` is NA, returns an sf object (filtered or unfiltered
#' based on `filter_return`). If `outfile` is a character string, returns
#' the character string of the saved file path.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#'  # load chm
#'  chm_fnm <- system.file(package = "cloud2trees", "extdata", "piles_chm.tif")
#'  my_chm <- terra::rast(chm_fnm)
#'  my_chm
#'  # load rgb
#'  rgb_fnm <- system.file(package = "cloud2trees", "extdata", "piles_rgb.tif")
#'  my_rgb <- terra::rast(rgb_fnm)
#'  my_rgb
#'  # where to write output?
#'  my_fnm <- tempfile()
#'  # piles_workflow()
#'  piles_workflow_ans <- piles_workflow(
#'    chm_rast = my_chm
#'    , seg_method = "dbscan"
#'    , min_ht_m = 1
#'    , max_ht_m = 6.6
#'    , min_area_m2 = pi*(1.5/2)^2
#'    , max_area_m2 = pi*(8/2)^2
#'    , min_convexity_ratio = 0.11
#'    , min_circularity_ratio = 0.22
#'    , smooth_segs = T
#'    , outfile = my_fnm
#'    , rgb_rast = my_rgb
#'    , red_band_idx = 1
#'    , green_band_idx = 2
#'    , blue_band_idx = 3
#'    , spectral_weight = 5
#'    , filter_return = F # don't filter the return so we can see what did not pass
#'  )
#'  # what does it have?
#'  names(piles_workflow_ans)
#'  # pile polygons
#'  piles_workflow_ans$segs_sf %>% dplyr::glimpse()
#'  # segmentation method parameters used
#'  piles_workflow_ans$seg_mthd_params
#'  # chm + structural piles since we didn't filter return
#'  terra::plot(my_chm, axes = F)
#'  piles_workflow_ans$segs_sf %>%
#'    terra::vect() %>%
#'    terra::plot(border = "magenta", lwd = 2, col = NA, add = T)
#'  # spectral threshold voting
#'  piles_workflow_ans$segs_sf %>%
#'    sf::st_drop_geometry() %>%
#'    dplyr::count(inrange_th_votes)
#'  # which would be removed if we set:
#'  # spectral_weight = 5 ?
#'  terra::plotRGB(
#'    my_rgb
#'    , main = "candidate piles kept (magenta) or removed (red)"
#'    , mar = c(0.2,0.2,2,0.2)
#'  )
#'  piles_workflow_ans$segs_sf %>%
#'    dplyr::filter(inrange_th_votes>=5) %>%
#'    terra::vect() %>%
#'    terra::plot(border = "magenta", lwd = 2, col = NA, add = T)
#'  piles_workflow_ans$segs_sf %>%
#'    dplyr::filter(inrange_th_votes<5) %>%
#'    terra::vect() %>%
#'    terra::plot(border = "red", lwd = 2, col = NA, add = T)
#' }
#'
piles_workflow <- function(
  # # piles_detect
  chm_rast
  , seg_method = "dbscan"
  , min_ht_m # set the min expected pile height
  , max_ht_m # set the max expected pile height
  , min_area_m2 # set the min expected pile area
  , max_area_m2 # set the max expected pile area
  , min_convexity_ratio # min required overlap between the predicted pile and the convex hull of the predicted pile
  , min_circularity_ratio
  , smooth_segs = T
  # # piles_spectral_filter
  , rgb_rast
  , red_band_idx
  , green_band_idx
  , blue_band_idx
  , spectral_weight = 4
  , filter_return = T
  , outfile = NA
) {
  check_segmentation_method_ans <- check_segmentation_method(method = seg_method)
  if(length(check_segmentation_method_ans)>1){
    warning(
      paste0( "...using only ", check_segmentation_method_ans[1], " method for segmentation")
    )
  }
  check_segmentation_method_ans <- check_segmentation_method_ans[1]
  ##############################################
  # piles_spectral_filter
  # checks before we get into it
  ##############################################
    if(!inherits(filter_return,"logical")){
      stop("Input `filter_return` should be logical: T to filter return based on the `spectral_weight` or F to return the full structural candidates")
    }
    chk_rgb_temp <- check_rgb_raster_bands(
      rast = rgb_rast
      , red_band_idx = red_band_idx
      , green_band_idx = green_band_idx
      , blue_band_idx = blue_band_idx
    )
    spectral_weight <- as.numeric(spectral_weight)[1]
    if(
      filter_return &&
      (
        is.na(spectral_weight) ||
        is.null(spectral_weight) ||
        is.nan(spectral_weight) ||
        !(spectral_weight %in% c(0:6))
      )
    ){
      stop("Input `spectral_weight` must be a number between 0 (no filtering based on spectral) and 6 (highest weighting of spectral data)")
    }
  ##############################################
  # piles_detect
  # checks will happen during
  ##############################################
    piles_detect_ans <- piles_detect(
      chm_rast = chm_rast
      , seg_method = check_segmentation_method_ans
      , min_ht_m = min_ht_m
      , max_ht_m = max_ht_m
      , min_area_m2 = min_area_m2
      , max_area_m2 = max_area_m2
      , min_convexity_ratio = min_convexity_ratio
      , min_circularity_ratio = min_circularity_ratio
      , smooth_segs = smooth_segs
      , outfile = outfile
    )
  ##############################################
  # write them if didn't already
  ##############################################
    if(
      inherits(piles_detect_ans,"list")
      && is.na(outfile)
    ){
      # get the things
      segs_sf <- piles_detect_ans[["segs_sf"]]
      seg_mthd_params <- piles_detect_ans[["seg_mthd_params"]]
      slice_chm_rast <- piles_detect_ans[["slice_chm_rast"]]
      # write the segs
      outfile <- "piles_workflow_piles.gpkg"
      sf::st_write(segs_sf, dsn = outfile, quiet = T, append = F)
      message(paste0("wrote structural candidate pile polygons to: \n      ",normalizePath(outfile)))
    }else if(file.exists(piles_detect_ans)){ # wrote the file in piles_detect
      # get the things
      outfile <- piles_detect_ans
      segs_sf <- sf::st_read(piles_detect_ans, quiet = T)
      seg_mthd_params <- get_segmentation_params(
        max_ht_m = max_ht_m
        , min_ht_m = min_ht_m
        , min_area_m2 = min_area_m2
        , max_area_m2 = max_area_m2
        , rast_res_m = terra::res(chm_rast)[1]
      )
      slice_chm_rast <- NA
    }else{
      stop("piles_detect() did not identify piles or could not write to provided outfile")
    }
  ##############################################
  # piles_spectral_filter
  ##############################################
    piles_spectral_filter_ans <- piles_spectral_filter(
      sf_data = segs_sf
      , rgb_rast = rgb_rast
      , red_band_idx = red_band_idx
      , green_band_idx = green_band_idx
      , blue_band_idx = blue_band_idx
      , spectral_weight = spectral_weight
      , filter_return = filter_return
    )
  # write
    sf::st_write(piles_spectral_filter_ans[["segs_sf"]], dsn = outfile, quiet = T, append = F)
    message(paste0("wrote spectrally filtered pile polygons to: \n      ",normalizePath(outfile)))
  # return
    return(list(
      segs_sf = piles_spectral_filter_ans[["segs_sf"]]
      , seg_mthd_params = seg_mthd_params
      , slice_chm_rast = slice_chm_rast
      , rgb_indices_rast = piles_spectral_filter_ans[["rgb_indices_rast"]]
      , outfile = normalizePath(outfile)
    ))
}
