testthat::test_that("itd_tuning() returns a list of objects", {
  ## test ws_fn_list
  testthat::expect_type(
    object = suppressWarnings({
      itd_tuning(
        input_las_dir = system.file(package = "lidR", "extdata", "MixedConifer.laz")
        , n_samples = 2
      ) %>% purrr::pluck("ws_fn_list")
    })
    , type = "list"
  )
  ## test plot_samples
  testthat::expect_s3_class(
    object = suppressWarnings({
      itd_tuning(
        input_las_dir = system.file(package = "lidR", "extdata", "MixedConifer.laz")
        , n_samples = 2
      ) %>% purrr::pluck("plot_samples")
    })
    , class = c("ggplot")
  )
})
