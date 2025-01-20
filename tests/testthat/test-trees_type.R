testthat::test_that("trees_type() returns a list of objects", {
  ## test class
  testthat::expect_type(
    object = suppressWarnings({
      trees_type(
        tree_list = dplyr::tibble(
          treeID = c(1:21)
          , tree_x = rnorm(n=21, mean = 458064, sd = 11)
          , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
        )
        , crs = "32613"
      )
    })
    , type = "list"
  )
  ## test tree_list
  testthat::expect_s3_class(
    object = suppressWarnings({
      trees_type(
        tree_list = dplyr::tibble(
          treeID = c(1:21)
          , tree_x = rnorm(n=21, mean = 458064, sd = 11)
          , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
        )
        , crs = "32613"
      ) %>% purrr::pluck("tree_list")
    })
    , class = c("sf", "data.frame")
  )
  ## test tree_list
  testthat::expect_s4_class(
    object = suppressWarnings({
      trees_type(
        tree_list = dplyr::tibble(
          treeID = c(1:21)
          , tree_x = rnorm(n=21, mean = 458064, sd = 11)
          , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
        )
        , crs = "32613"
      ) %>% purrr::pluck("foresttype_rast")
    })
    , class = c("SpatRaster")
  )
  ## test error message
  testthat::expect_error(
    object = suppressWarnings({
      trees_type(
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
      trees_type(
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
