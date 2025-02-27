testthat::test_that("itd_ws_functions() returns a list of functions", {
  ## test class
  testthat::expect_type(
    object = suppressWarnings({
      itd_ws_functions()
    })
    , type = "list"
  )
  ## test value
  testthat::expect_true(
    object = suppressWarnings({
      itd_ws_functions()[["exp_fn"]] %>% is.function()
    })
  )

})
