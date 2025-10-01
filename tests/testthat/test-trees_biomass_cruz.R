testthat::test_that("trees_biomass_cruz() returns a list of objects", {
  suppressWarnings(suppressMessages({
      trees_biomass_cruz_ans <- trees_biomass_cruz(
        tree_list = dplyr::tibble(
          treeID = c(1:21)
          , tree_x = rnorm(n=21, mean = 458064, sd = 11)
          , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
          , tree_height_m = rep(6,times=21)
          , crown_area_m2 = rep(2,times=21)
          , dbh_cm = rep(10,times=21)
          , tree_cbh_m = rep(2,times=21)
        )
        , crs = "32613"
      )
    }))
  ####################
  # tests
  ####################
  ## test class
  testthat::expect_type(
    object = trees_biomass_cruz_ans
    , type = "list"
  )

  ## includes the named list?
  testthat::expect_named(trees_biomass_cruz_ans, expected = c("tree_list", "stand_cell_data"), ignore.order = TRUE)

  ## test class
  testthat::expect_type(
    object = trees_biomass_cruz_ans
    , type = "list"
  )
  ## test tree_list
  testthat::expect_s3_class(
    object = trees_biomass_cruz_ans %>% purrr::pluck("tree_list")
    , class = c("sf", "data.frame")
  )
  ## test stand_cell_data
  testthat::expect_s3_class(
    object = trees_biomass_cruz_ans %>% purrr::pluck("stand_cell_data")
    , class = c("sf", "data.frame")
  )
  ## test error message
  testthat::expect_error(
    object = suppressWarnings({
      trees_biomass_cruz(
        tree_list = dplyr::tibble(
          treeID = c(1:21)
          , tree_x = rnorm(n=21, mean = 458064, sd = 11)
          , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
          , tree_height_m = rep(6,times=21)
          , crown_area_m2 = rep(2,times=21)
          # , dbh_cm = rep(10,times=21)
          # , tree_cbh_m = rep(2,times=21)
        )
        , crs = "32613"
      )
    })
  )

})
