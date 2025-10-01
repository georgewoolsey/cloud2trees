testthat::test_that("trees_type() returns a list of objects", {
  suppressWarnings(suppressMessages({
    trees_type_ans <- trees_type(
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
  ## includes the named list?
  testthat::expect_named(trees_type_ans, expected = c("tree_list", "foresttype_rast"), ignore.order = TRUE)

  ## test class
  testthat::expect_type(
    trees_type_ans
    , type = "list"
  )
  ## test tree_list
  testthat::expect_s3_class(
    object = trees_type_ans %>% purrr::pluck("tree_list")
    , class = c("sf", "data.frame")
  )
  ## test tree_list
  testthat::expect_s4_class(
    object = trees_type_ans %>% purrr::pluck("foresttype_rast")
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
