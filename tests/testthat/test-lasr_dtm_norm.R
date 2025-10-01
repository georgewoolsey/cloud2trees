testthat::test_that("lasr_dtm_norm() returns a LASRpipeline", {
  testthat::expect_s3_class(
    lasr_dtm_norm("testtest.tif")
    # , c("LASRpipeline", "list") # lasR 0.13.1 and below
    , c("PipelinePtr") # lasR at least 0.17.2 maybe earlier
  )
})
