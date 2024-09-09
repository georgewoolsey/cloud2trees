testthat::test_that("chunk_las_catalog() returns a list with 4 named items", {
  testthat::expect_length(
    object = chunk_las_catalog(
      folder = list.files(system.file(package = "lidR", "extdata/"), recursive = T, full.names = T) %>%
        tolower() %>%
        stringr::str_subset("mixedconifer") %>%
        stringr::str_subset(".*\\.(laz|las)$")
      , outfolder = getwd()
    )
    , n = 4
  )
  testthat::expect_named(
    object = chunk_las_catalog(
      folder = list.files(system.file(package = "lidR", "extdata/"), recursive = T, full.names = T) %>%
        tolower() %>%
        stringr::str_subset("mixedconifer") %>%
        stringr::str_subset(".*\\.(laz|las)$")
      , outfolder = getwd()
    )
    , expected = c("process_data", "grid_subset_switch", "plt", "las_ctg")
  )
})
