testthat::test_that("cloud2raster() returns a list with 5 named items", {
  testthat::expect_named(
    object = cloud2raster(
      input_las_dir = list.files(system.file(package = "lidR", "extdata/"), recursive = T, full.names = T) %>%
        tolower() %>%
        stringr::str_subset("mixedconifer") %>%
        stringr::str_subset(".*\\.(laz|las)$")
      , output_dir = tempdir()
    )
    , expected = c("dtm_rast"
      , "chm_rast"
      , "create_project_structure_ans"
      , "chunk_las_catalog_ans"
      , "normalize_flist")
  )
})
