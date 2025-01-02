#' use lasR to combine DTM and normalize step
#' @description
#' Combining DTM and normalize step using `lasR` functionality
#'
#' @param dtm_file_name string. Where to write the DTM.
#' @param frac_for_tri numeric. The fraction of points used in Delauny triangulation.
#' @param dtm_res numeric. The desired resolution of the DTM produced in meters.
#' @param norm_accuracy numeric. see `chunk_las_catalog`. Choose processing accuracy.
#'      accuracy_level = 1 uses DTM to height normalize the points
#'      accuracy_level = 2 uses triangulation with high point density (20 pts/m2) to height normalize the points
#'      accuracy_level = 3 uses triangulation with very high point density (100 pts/m2) to height normalize the points
#'
#' @return A `lasR` pipeline
#'
#' @references
#' [https://r-lidar.github.io/lidRbook/normalization.html](https://r-lidar.github.io/lidRbook/normalization.html)
#' [https://github.com/r-lidar/lasR/issues/18#issuecomment-2027818414](https://github.com/r-lidar/lasR/issues/18#issuecomment-2027818414)
#'
#' @keywords internal
#'
lasr_dtm_norm <- function(
  dtm_file_name
  , frac_for_tri = 1
  , dtm_res = 1
  , norm_accuracy = 2
){
  # lasR required
  if(!requireNamespace("lasR", quietly = TRUE)) {
    stop(paste0(
      "Package \"lasR\" must be installed to use this function."
      , "\n"
      , "try `install.packages(\"lasR\", repos = \"https://r-lidar.r-universe.dev\")`"
    ))
  }
  # perform Delaunay triangulation
    # tri = lasR::triangulate(filter = "-keep_class 2 -keep_class 9 -keep_random_fraction 0.01")
    ####
    # set filter based on # points
    ####
    ################################
    ## lasR 0.13 update to decimate ground points for triangulation
    ## 2024-03-29 version: https://github.com/r-lidar/lasR/issues/18
    ## 2024-12-18 version:: https://github.com/r-lidar/lasR/issues/102
    ################################
    # resample ground for super high density clouds
    if( dplyr::coalesce(as.numeric(frac_for_tri),1)<1 ){
      samp_res <- dplyr::case_when(
        as.numeric(frac_for_tri) <= 0.01 ~ 0.4 # 40 cm
        , as.numeric(frac_for_tri) <= 0.3 ~ 0.2 # 20 cm
        , T ~ 0.1 # 10 cm
      )
      lasr_resample_gnd <- lasR::sampling_pixel(
        res = samp_res
        , method = "min"
        , filter = lasR::keep_ground_and_water()
      )
    }else{ # fake the stage
      lasr_resample_gnd <- lasR::summarise()
    }
    filter_for_dtm <- lasR::keep_ground_and_water()

    ####
    # triangulate with filter
    # produces a triangulation of the ground points (meshed DTM)
    ####
    lasr_triangulate <- lasR::triangulate(
      # class 2 = ground; class 9 = water
      filter = filter_for_dtm
      , max_edge = 0
      # , max_edge = c(0,1)
      # # write to disk to preserve memory
      , ofile = ""
      # , ofile = paste0(config$las_denoise_dir, "/", "*_tri.gpkg")
    )
  # rasterize the result of the Delaunay triangulation
    lasr_dtm <- lasR::rasterize(
      res = dtm_res
      , operators = lasr_triangulate
      , filter = lasR::drop_noise()
      # # write to disk to preserve memory
      , ofile = dtm_file_name
    )
  # normalize
  # stage = triangulate: takes foreevvveerrrrrrr
    # ... but see: https://github.com/r-lidar/lasR/issues/18#issuecomment-2027818414
    ## at this density of point, my advice is anyway to decimate your ground points.
    ## With 1/100 of the ground points you already have 5 ground pts/m2 to compute a
    ## DTM with a 1 m resolution! You could even decimate to 1/250.
    ## This will solve your computation time issue in the same time.
    ## stage = dtm
      # ... see: https://github.com/r-lidar/lasR/issues/17#issuecomment-2027698100
      # also from the lidR book https://r-lidar.github.io/lidRbook/norm.html:
        ## "Point cloud normalization without a DTM interpolates the elevation of
        ## every single point locations using ground points. It no longer uses elevations
        ## at discrete predefined locations. Thus the methods is exact, computationally speaking.
        ## It means that it is equivalent to using a continuous DTM but it is important
        ## to recall that all interpolation methods are interpolation and by definition
        ## make guesses with different strategies. Thus by “exact” we mean “continuous”.
    if(as.numeric(norm_accuracy) %in% c(2,3)){
      lasr_normalize <- lasR::transform_with(
        stage = lasr_triangulate
        , operator = "-"
      )
    }else{
      lasr_normalize <- lasR::transform_with(
        stage = lasr_dtm
        , operator = "-"
      )
    }
  # pipeline
  pipeline <- lasr_resample_gnd + lasr_triangulate + lasr_dtm + lasr_normalize
  return(pipeline)
}
