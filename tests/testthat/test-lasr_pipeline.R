testthat::test_that("lasr_pipeline() returns a list", {
  testthat::expect_type(
    object = lasr_pipeline(
      processing_grid_num = 1
      , process_data = chunk_las_catalog(
          folder = system.file("extdata", "MixedConifer.laz", package="lidR")
          , outfolder = tempdir()
        )[["process_data"]]
    )
    , type = "list"
  )
})
