testthat::test_that("simplify_multipolygon_crowns() returns a spatial data frame", {
  testthat::expect_s3_class(
    object = simplify_multipolygon_crowns(
      sf::st_read( paste0(system.file(package = "cloud2trees"),"/extdata/crowns_poly.gpkg"), quiet = T)
    )
    , class = c("sf","data.frame")
  )
})
