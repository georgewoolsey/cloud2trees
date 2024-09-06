testthat::test_that("lasr_dtm_norm() returns a LASRpipeline", {
  testthat::expect_s3_class(
    lasr_dtm_norm("testtest.tif")
    , c("LASRpipeline", "list")
  )
})
