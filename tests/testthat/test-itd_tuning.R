testthat::test_that("itd_tuning() returns a list of objects", {
  suppressWarnings(suppressMessages({
      itd_tuning_ans <- itd_tuning(
        input_las_dir = system.file(package = "lidR", "extdata", "MixedConifer.laz")
        , n_samples = 2
      )
    }))
  ####################
  # tests
  ####################
  ## includes the named list?
  testthat::expect_named(itd_tuning_ans, expected = c("ws_fn_list", "plot_samples", "plot_sample_summary", "crowns"), ignore.order = TRUE)
  # get the things
  ws_fn_list <- itd_tuning_ans %>% purrr::pluck("ws_fn_list")
  plot_samples <- itd_tuning_ans %>% purrr::pluck("plot_samples")
  plot_sample_summary <- itd_tuning_ans %>% purrr::pluck("plot_sample_summary")
  ## test ws_fn_list
  testthat::expect_type(
    object = ws_fn_list
    , type = "list"
  )
  ## test plot_samples
  testthat::expect_s3_class(
    object = plot_samples
    , class = c("ggplot")
  )
  ## test plot_sample_summary
  testthat::expect_s3_class(
    object = plot_sample_summary
    , class = c("ggplot")
  )

  ####################
  # CHM
  ####################
  suppressWarnings(suppressMessages({
      itd_tuning_ans2 <- itd_tuning(
        input_chm_rast = terra::rast( system.file(package = "cloud2trees", "extdata", "chm.tif") )
        , n_samples = 2
      )
    }))
  ####################
  # tests
  ####################
  ## includes the named list?
  testthat::expect_named(itd_tuning_ans2, expected = c("ws_fn_list", "plot_samples", "plot_sample_summary", "crowns"), ignore.order = TRUE)
  # get the things
  ws_fn_list2 <- itd_tuning_ans2 %>% purrr::pluck("ws_fn_list")
  plot_samples2 <- itd_tuning_ans2 %>% purrr::pluck("plot_samples")
  plot_sample_summary2 <- itd_tuning_ans2 %>% purrr::pluck("plot_sample_summary")
  ## test ws_fn_list
  testthat::expect_type(
    object = ws_fn_list2
    , type = "list"
  )
  ## test plot_samples
  testthat::expect_s3_class(
    object = plot_samples2
    , class = c("ggplot")
  )
  ## test plot_sample_summary
  testthat::expect_s3_class(
    object = plot_sample_summary2
    , class = c("ggplot")
  )
})
