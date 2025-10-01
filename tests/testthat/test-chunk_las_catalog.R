testthat::test_that("chunk_las_catalog() returns a list with 4 named items", {
  suppressWarnings(suppressMessages({
      chunk_las_catalog_ans <- chunk_las_catalog(
        folder = system.file("extdata", "MixedConifer.laz", package="lidR")
        , outfolder = tempdir()
      )
    }))
  ####################
  # tests
  ####################
  ## test class
  testthat::expect_type(
    object = chunk_las_catalog_ans
    , type = "list"
  )

  ## includes the named list?
  testthat::expect_named(chunk_las_catalog_ans, expected = c("process_data", "is_chunked_grid", "plt", "las_ctg"), ignore.order = TRUE)
  ## length
  testthat::expect_length(
    object = chunk_las_catalog_ans
    , n = 4
  )

  ## test process_data
  testthat::expect_s3_class(
    object = chunk_las_catalog_ans %>% purrr::pluck("process_data")
    , class = c("sf", "data.frame")
  )

  ## test plt
  testthat::expect_s3_class(
    object = chunk_las_catalog_ans$plt
    , class = c("ggplot")
  )

})
