testthat::test_that("lasr_chm() returns a LASRpipeline", {
  testthat::expect_s3_class(
    lasr_chm(chm_file_name = "testtest.tif")
    , c("LASRpipeline", "list")
  )
})
