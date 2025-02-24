testthat::test_that("rf_tune_subsample() returns a number for mtry", {
  ## test class
  testthat::expect_type(
    object = suppressWarnings({
      rf_tune_subsample(
        predictors = datasets::trees[,-1]
        , response = datasets::trees[,1]
      ) # %>% typeof()
    })
    , type = "double"
  )
  ## test value
  testthat::expect_length(
    object = suppressWarnings({
      rf_tune_subsample(
        predictors = datasets::trees[,-1]
        , response = datasets::trees[,1]
      )
    })
    , 1
  )

})
