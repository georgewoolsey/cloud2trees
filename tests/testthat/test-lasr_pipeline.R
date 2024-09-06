testthat::test_that("lasr_pipeline() returns a list", {
  testthat::expect_type(
    object = lasr_pipeline(
      processing_grid_num = 1
      , process_data = chunk_las_catalog(
          folder =
            # "C:/data/usfs/point_cloud_tree_detection_ex/data/testtest"
            list.files(system.file(package = "lidR", "extdata/"), recursive = T, full.names = T) %>%
              tolower() %>%
              stringr::str_subset("mixedconifer")
          , outfolder = getwd()
        )[["process_data"]]
    )
    , type = "list"
  )
})
