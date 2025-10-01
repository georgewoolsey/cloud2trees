testthat::test_that("trees_landfire_cbd() returns a list of objects", {
  suppressWarnings(suppressMessages({
      trees_landfire_cbd_ans <- trees_landfire_cbd(
        tree_list = dplyr::tibble(
          treeID = c(1:21)
          , tree_x = rnorm(n=21, mean = 458064, sd = 11)
          , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
        )
        , crs = "32613"
      )
    }))
  ####################
  # tests
  ####################
  ## test class
  testthat::expect_type(
    object = trees_landfire_cbd_ans
    , type = "list"
  )

  ## includes the named list?
  testthat::expect_named(trees_landfire_cbd_ans, expected = c("tree_list", "landfire_rast"), ignore.order = TRUE)

  ## test tree_list
  testthat::expect_s3_class(
    object = trees_landfire_cbd_ans %>% purrr::pluck("tree_list")
    , class = c("sf", "data.frame")
  )
  ## test tree_list
  testthat::expect_s4_class(
    object = trees_landfire_cbd_ans %>% purrr::pluck("landfire_rast")
    , class = c("SpatRaster")
  )
  ## test error message
  testthat::expect_error(
    object = suppressWarnings({
      trees_landfire_cbd(
        tree_list = dplyr::tibble(
          treeID = c(1:21)
          , tree_x = rnorm(n=21, mean = 458064, sd = 11)
          , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
        )
        # , crs = "32613" # took out crs
      )
    })
  )
  ## test no error message
  testthat::expect_no_error(
    object = suppressWarnings({
      trees_landfire_cbd(
        tree_list = dplyr::tibble(
          treeID = c(1:21)
          , tree_x = rnorm(n=21, mean = 458064, sd = 11)
          , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
        )
        , crs = "32613"
      )
    })
  )

})
