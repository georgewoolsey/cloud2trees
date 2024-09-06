testthat::test_that("chunk_las_catalog() returns a list with three named items", {
  testthat::expect_length(
    object = chunk_las_catalog(
      folder = list.files(system.file(package = "lidR", "extdata/"), recursive = T, full.names = T) %>%
        tolower() %>%
        stringr::str_subset("mixedconifer")
      , outfolder = getwd()
    )
    , n = 3
  )
  testthat::expect_named(
    object = chunk_las_catalog(
      folder = list.files(system.file(package = "lidR", "extdata/"), recursive = T, full.names = T) %>%
        tolower() %>%
        stringr::str_subset("mixedconifer")
      , outfolder = getwd()
    )
    , expected = c("process_data", "grid_subset_switch", "plt")
  )
})
