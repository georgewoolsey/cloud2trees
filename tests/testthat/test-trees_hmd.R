testthat::test_that("trees_hmd() returns a spatial data frame", {
  testthat::expect_s3_class(
    object = suppressWarnings({
      trees_hmd(
        trees_poly = sf::st_read( paste0(system.file(package = "cloud2trees"),"/extdata/crowns_poly.gpkg"), quiet = T)
        , norm_las = paste0(system.file(package = "cloud2trees"),"/extdata/norm_las")
        , estimate_missing_hmd = T
        , force_same_crs = T
      )
    })
    , class = c("sf", "data.frame")
  )

})
