testthat::test_that("treels_stem_dbh() returns a data.frame", {
  testthat::expect_s3_class(
    object = treels_stem_dbh(
      folder = system.file("extdata", "MixedConifer.laz", package="lidR")
      , outfolder = getwd()
    )
    , class = c("sf", "data.frame")
  )
})
