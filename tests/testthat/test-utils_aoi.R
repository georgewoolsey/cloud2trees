testthat::test_that("functions in utils_aoi handle data correctly for cloud2trees_to_lanl_trees()", {
  # load data
  suppressWarnings(suppressMessages({
    f <- system.file("extdata", "crowns_poly.gpkg", package = "cloud2trees")
    tree_list <- sf::st_read(f,quiet=T)
    aoi <- tree_list %>% dplyr::slice_sample(n=1) %>% sf::st_centroid() %>% sf::st_buffer(20)
    # ggplot2::ggplot() + ggplot2::geom_sf(data=aoi) + ggplot2::geom_sf(data=tree_list)
  }))
  ###################################################
  # perform tests on the processed_data
  ###################################################

  #######################################
  # clip_tree_list_aoi
  #######################################
    suppressWarnings(suppressMessages({
      clip_tree_list_aoi_ans <- clip_tree_list_aoi(
        tree_list = tree_list
        , study_boundary = aoi
        , bbox_aoi = F
        , buffer = 0
        , reproject_epsg = NULL
      )
    }))
    ####################
    # tests
    ####################
    ## test class
    testthat::expect_type(object = clip_tree_list_aoi_ans, type = "list")
    ## test objects
    testthat::expect_s3_class(object = clip_tree_list_aoi_ans$tree_list, class = c("sf", "data.frame"))
    testthat::expect_s3_class(object = clip_tree_list_aoi_ans$aoi, class = c("sf", "data.frame"))
    testthat::expect_lt(object = nrow(clip_tree_list_aoi_ans$tree_list), expected = nrow(tree_list))

  #######################################
  # quicfire_define_domain
  #######################################
    suppressWarnings(suppressMessages({
      quicfire_domain_df <- quicfire_define_domain(
        sf_data = clip_tree_list_aoi_ans$aoi
        , horizontal_resolution = 2
      )
    }))

    ####################
    # tests
    ####################
    ## test class
    testthat::expect_type(object = quicfire_domain_df, type = "list")
    ## test objects
    testthat::expect_s3_class(object = quicfire_domain_df$quicfire_domain_df, class = c("sf", "data.frame"))
    testthat::expect_type(object = quicfire_domain_df$domain_path, type = "character")
    testthat::expect_equal(object = nrow(quicfire_domain_df$quicfire_domain_df), expected = 1)

  #######################################
  # quicfire_dtm_topofile
  #######################################
    dtmf <- system.file("extdata", "chm.tif", package = "cloud2trees") # fake a dtm...just need a rast
    # dtm <- terra::rast(dtmf)
      # terra::plot(dtm)
      # terra::plot(
      #   tree_list %>%
      #     sf::st_transform(terra::crs(dtm)) %>%
      #     terra::vect()
      #   , add = T, border = "red", col = NA
      # )

    suppressWarnings(suppressMessages({
      quicfire_dtm_topofile_ans <- quicfire_dtm_topofile(
        dtm_rast = dtmf
        , study_boundary = quicfire_domain_df$quicfire_domain_df
      )
    }))

    ####################
    # tests
    ####################
    ## test class
    testthat::expect_type(object = quicfire_dtm_topofile_ans, type = "list")
    ## test objects
    testthat::expect_s4_class(object = quicfire_dtm_topofile_ans$dtm, class = "SpatRaster")

  #######################################
  # make_lanl_trees_input
  #######################################
    suppressWarnings(suppressMessages({
      make_lanl_trees_input_ans <- make_lanl_trees_input(
        tree_list = clip_tree_list_aoi_ans$tree_list %>%
          dplyr::mutate(crown_dia_m = 21, max_crown_diam_height_m = 3, landfire_tree_kg_per_m3=4)
        , quicfire_domain_df = quicfire_domain_df$quicfire_domain_df
        , topofile = "flat"
        , cbd_col_name = "landfire_tree_kg_per_m3"
      )
    }))

    ####################
    # tests
    ####################
    ## test class
    testthat::expect_type(object = make_lanl_trees_input_ans, type = "list")
    ## test objects
    testthat::expect_s3_class(object = make_lanl_trees_input_ans$treelist, class = "data.frame")
    testthat::expect_type(object = make_lanl_trees_input_ans$fuellist_path, type = "character")
    testthat::expect_type(object = make_lanl_trees_input_ans$treelist_path, type = "character")

    testthat::expect_error(
      object =
        make_lanl_trees_input(
          tree_list = clip_tree_list_aoi_ans$tree_list
          , quicfire_domain_df = quicfire_domain_df$quicfire_domain_df
          , topofile = "flat"
          , cbd_col_name = "landfire_tree_kg_per_m3"
        )
      # , "landfire_tree_kg_per_m3 doesn't even exist guy"
    )


})

