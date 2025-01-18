testthat::test_that("trees_biomass_cruz() returns a list of objects", {
  ## test class
  testthat::expect_type(
    object = suppressWarnings({
      trees_biomass_cruz(
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
    })
    , type = "list"
  )
  ## test tree_list
  testthat::expect_s3_class(
    object = suppressWarnings({
      trees_biomass_cruz(
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
      ) %>% purrr::pluck("tree_list")
    })
    , class = c("sf", "data.frame")
  )
  ## test stand_cell_data
  testthat::expect_s3_class(
    object = suppressWarnings({
      trees_biomass_cruz(
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
      ) %>% purrr::pluck("stand_cell_data")
    })
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
  ## test no error message
  testthat::expect_no_error(
    object = suppressWarnings({
      trees_biomass_cruz(
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
    })
  )

})
