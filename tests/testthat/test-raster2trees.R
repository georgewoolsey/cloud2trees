testthat::test_that("raster2trees() returns a spatial data frame", {
  testthat::expect_s3_class(
    object = raster2trees(
      chm_rast = terra::rast( paste0(system.file(package = "cloud2trees"),"/extdata/chm.tif") )
      , outfolder = tempdir()
    )
    , class = c("sf","data.frame")
  )
})
