testthat::test_that("trees_dbh() returns a data.frame", {
  ## test error message
  testthat::expect_error(
    object = suppressWarnings({
      trees_dbh(
        tree_list = dplyr::tibble(
          treeID = c(1:21)
          , tree_x = rnorm(n=21, mean = 458064, sd = 11)
          , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
          , tree_height_m = exp(rgamma(n = 21, shape = (7/4)^2, rate = (4^2)/7))
        )
        ## , crs = "32613" # took out crs
      )
    })
  )
  ## test error message
  testthat::expect_error(
    object = suppressWarnings({
      trees_dbh(
        tree_list = dplyr::tibble(
          treeID = c(1:21)
          , tree_x = rnorm(n=21, mean = 458064, sd = 11)
          , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
          ## took out tree_height_m
        )
        , crs = "32613"
      )
    })
  )
})
