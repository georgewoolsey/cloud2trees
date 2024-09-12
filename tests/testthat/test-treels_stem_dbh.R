testthat::test_that("treels_stem_dbh() returns a data.frame", {
  testthat::expect_s3_class(
    object = treels_stem_dbh(
      folder = list.files(system.file(package = "lidR", "extdata/"), recursive = T, full.names = T) %>%
        tolower() %>%
        stringr::str_subset("mixedconifer") %>%
        stringr::str_subset(".*\\.(laz|las)$")
      , outfolder = getwd()
    )
    , class = "data.frame"
  )
})
