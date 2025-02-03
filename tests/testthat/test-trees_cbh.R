testthat::test_that("trees_cbh() returns a spatial data frame", {
  testthat::expect_s3_class(
    object = trees_cbh(
      trees_poly = sf::st_read( paste0(system.file(package = "cloud2trees"),"/extdata/crowns_poly.gpkg"), quiet = T)
      , norm_las = paste0(system.file(package = "cloud2trees"),"/extdata/norm_las")
      , tree_sample_n = 30
      , estimate_missing_cbh = T
      , force_same_crs = T
    )
    , class = c("sf","data.frame")
  )
})
