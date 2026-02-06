testthat::test_that("piles_detect() returns a list of objects. piles_spectral_filter() does work.", {
  ###########################################
  # piles_detect
  ###########################################
  suppressWarnings(suppressMessages({
      my_chm <- terra::rast(system.file(package = "cloud2trees", "extdata", "piles_chm.tif"))
      piles_detect_ans <- piles_detect(
        chm_rast = my_chm
        , seg_method = "dbscan"
        , min_ht_m = 1
        , max_ht_m = 4
        , min_area_m2 = pi*(1.5/2)^2
        , max_area_m2 = pi*(6/2)^2
        , min_convexity_ratio = 0.3
        , min_circularity_ratio = 0.4
        , smooth_segs = T
        , outfile = NA
      )
    }))
  ####################
  # tests
  ####################
  ## includes the named list?
  testthat::expect_named(piles_detect_ans, expected = c("segs_sf", "seg_mthd_params", "slice_chm_rast"), ignore.order = TRUE)
  # get the things
  segs_sf <- piles_detect_ans[["segs_sf"]]
  seg_mthd_params <- piles_detect_ans[["seg_mthd_params"]]
  slice_chm_rast <- piles_detect_ans[["slice_chm_rast"]]
  ## test segs_sf
  testthat::expect_s3_class(
    object = segs_sf
    , class = "sf"
  )
  testthat::expect_s3_class(
    object = segs_sf
    , class = "data.frame"
  )
  ## test seg_mthd_params
  testthat::expect_type(
    object = seg_mthd_params
    ,type = "list"
  )
  ## test slice_chm_rast
  testthat::expect_s4_class(
    object = slice_chm_rast
    , "SpatRaster"
  )
  ####################
  # others
  ####################
  ## test error
  testthat::expect_error(
    object = suppressWarnings({
      piles_detect(
        chm_rast = my_chm
        , seg_method = "ferp"
        , min_ht_m = 1
        , max_ht_m = 4
        , min_area_m2 = pi*(1.5/2)^2
        , max_area_m2 = pi*(6/2)^2
        , min_convexity_ratio = 0.3
        , min_circularity_ratio = 0.4
        , smooth_segs = T
        , outfile = NA
      )
    })
  )
  ## test success
  testthat::expect_no_error(
    suppressWarnings(suppressMessages({
      piles_detect(
        chm_rast = my_chm
        , seg_method = "watershed"
        , min_ht_m = 1
        , max_ht_m = 2.2
        , min_area_m2 = pi*(1.5/2)^2
        , max_area_m2 = pi*(4/2)^2
        , min_convexity_ratio = 0.3
        , min_circularity_ratio = 0.4
        , smooth_segs = T
        , outfile = NA
      )
    }))
  )
  ###########################################
  # piles_spectral_filter
  ###########################################
  suppressWarnings(suppressMessages({
      my_rgb <- terra::rast(system.file(package = "cloud2trees", "extdata", "piles_rgb.tif"))
      piles_spectral_filter_ans <- piles_spectral_filter(
        sf_data = segs_sf
        , rgb_rast = my_rgb
        , red_band_idx = 1
        , green_band_idx = 2
        , blue_band_idx = 3
        , spectral_weight = 5
        , filter_return = T
      )
    }))
  ####################
  # tests
  ####################
  ## includes the named list?
  testthat::expect_named(piles_spectral_filter_ans, expected = c("segs_sf", "rgb_indices_rast"), ignore.order = TRUE)
  # get the things
  segs_sf <- piles_spectral_filter_ans[["segs_sf"]]
  rgb_indices_rast <- piles_spectral_filter_ans[["rgb_indices_rast"]]
  ## test segs_sf
  testthat::expect_s3_class(
    object = segs_sf
    , class = "sf"
  )
  testthat::expect_s3_class(
    object = segs_sf
    , class = "data.frame"
  )
  ## test rgb_indices_rast
  testthat::expect_s4_class(
    object = rgb_indices_rast
    , "SpatRaster"
  )


})
