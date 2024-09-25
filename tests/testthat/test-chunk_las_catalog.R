testthat::test_that("chunk_las_catalog() returns a list with 4 named items", {
  testthat::expect_length(
    object = chunk_las_catalog(
      folder = system.file("extdata", "MixedConifer.laz", package="lidR")
      , outfolder = tempdir()
    )
    , n = 4
  )
  testthat::expect_named(
    object = chunk_las_catalog(
      folder = system.file("extdata", "MixedConifer.laz", package="lidR")
      , outfolder = tempdir()
    )
    , expected = c("process_data", "is_chunked_grid", "plt", "las_ctg")
  )
})
