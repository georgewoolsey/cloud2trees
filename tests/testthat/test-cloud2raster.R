testthat::test_that("cloud2raster() returns a list with 5 named items", {
  testthat::expect_named(
    object = cloud2raster(
      input_las_dir = system.file("extdata", "MixedConifer.laz", package="lidR")
      , output_dir = tempdir()
    )
    , expected = c("dtm_rast"
      , "chm_rast"
      , "create_project_structure_ans"
      , "chunk_las_catalog_ans"
      , "normalize_flist")
  )
})
