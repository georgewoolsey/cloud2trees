#' @title internal functions to extract raster values at point locations
#'
#' @description
#' internal functions to extract raster values at point locations used by [trees_type()] and tree_biomass() for LANDFIRE rasters
#' for the most part, these functions take point and raster data as inputs
#' the general process is (see [trees_type()] for example):
#'
#' * crop_raster_match_points()
#' * check for undesirable or NA values in the point_values from crop_raster_match_points()
#' * if undesirable values:
#'    - reclass_landfire_rast() or reclass_foresttype_rast()
#'    - agg_fill_rast_match_points(), see details above function definition
#'    - use the point_values from agg_fill_rast_match_points() if not null
#'    - otherwise, use the point_values from crop_raster_match_points()
#'
#' @param points sf.
#' @param rast SpatRaster.
#' @param study_boundary sf.
#' @param max_search_dist_m numeric.
#'
#' @keywords internal
#'
crop_raster_match_points <- function(
  points
  , rast
  , study_boundary = NA
  , max_search_dist_m = 1000
){
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! function 8
  # crop the raster to the points/study area + buffer
  # extract the raster values at the points
  # return the cropped raster and the points with the raster value
  # get extent of trees data
    bbox_temp <- points %>% dplyr::ungroup() %>% sf::st_bbox()
  # find largest side
    buffer_temp <- max(
        bbox_temp["xmax"]-bbox_temp["xmin"]
        , bbox_temp["ymax"]-bbox_temp["ymin"]
      ) %>%
      `*`(0.6) # reduce the buffer by 40%
      max(max_search_dist_m) # compare to user defined max

  # check against the study boundary
    if(inherits(study_boundary, "sf") || inherits(study_boundary, "sfc")){
      bbox_b_temp <- study_boundary %>%
        sf::st_union() %>%
        sf::st_as_sf() %>%
        sf::st_transform(sf::st_crs(tree_tops)) %>%
        sf::st_bbox()
      # find largest side and compare to current setting from trees
      buffer_b_temp <- max(
          bbox_b_temp["xmax"]-bbox_b_temp["xmin"]
          , bbox_b_temp["ymax"]-bbox_b_temp["ymin"]
        ) %>%
        `*`(0.6) # reduce the buffer by 40%
      # reset bbox and buffer if larger
      if(buffer_b_temp>buffer_temp){
        buffer_temp <- buffer_b_temp
        bbox_temp <- bbox_b_temp
      }
    }

  # apply the buffer to get the search extent
    ext_temp <- bbox_temp %>%
      sf::st_as_sfc() %>%
      sf::st_buffer(buffer_temp, endCapStyle = "SQUARE") %>%
      terra::vect() %>%
      terra::project(terra::crs(rast))

  ##################################
  # crop the raster to this extent
  ##################################
    rast_temp <- rast %>%
      terra::crop(ext_temp) # crop it first to make if faster

  ##################################
  # do we even need to fill NA values for this list?
  ##################################
    # attempt to extract the forest type for the trees
    # now we apply the `terra::extract()` function to attach forest type to our original tree list
    tree_ftype_temp <-
      terra::extract(
        x = rast_temp
        , y = points %>%
          terra::vect() %>%
          terra::project(terra::crs(rast)) # don't forget to reproject
        , ID = T # this way, the 2nd column will always be the raster value of the 1st layer
      ) %>%
      dplyr::rename(point_record_number = 1, raster_value = 2)
  # return
    return(list(
      point_values = tree_ftype_temp
      , rast = rast_temp
      , bbox = bbox_temp
    ))
}
#######################################################
# intermediate function 1
#######################################################
fill_rast_na <- function(rast){
    # fill NA values of raster using nearest neighbor interpolation
  # make sure there are filled cells
  if(
    terra::global(rast, fun = "isNA") == terra::ncell(rast)
  ){
    stop("All raster values are NA. Cannot fill entirely empty raster.")
  }

  # blank raster
  r_temp <- rast
  r_temp[] <- NA

  # get points with data
  p_temp <- terra::as.points(rast)

  # get nearest neighbor "prediction" as the values of the closest location with data using
  # terra::voronoi which creates a Voronoi diagram (also known as Delaunay triangles or Thiessen diagram)
  #   The Voronoi diagram is created when a region with n points is partitioned into convex polygons
  #   such that each polygon contains exactly one generating point, and every point in a given polygon
  #   is closer to its generating point than to any other.
  v_temp <- terra::voronoi(
    x = p_temp
    # set the boundary to the extent of the original raster
    , bnd = terra::ext(rast) %>% terra::vect()
  )

  # fill in our blank raster with these predictions
  r_temp <- terra::rasterize(x = v_temp, y = r_temp, field = names(v_temp)[1])

  # nearest neighbor interpolate applied via terra::cover
    # Replace NA or other values in SpatRaster x with the values of SpatRaster y
  rast_fill <- terra::cover(
    x = rast
    , y = r_temp
    , values = NA
  )

  # return
  return(rast_fill)
}
#######################################################
# intermediate function 2
#######################################################
  # determine factor to fill
  get_fact_fn <- function(a_m2,res=30) {
    if(res==30){
      # 150B = 2*3; 200B = 2*3;...;500B = 5*3; 550B = 6*3;...; 900B = 9*3
      fact <- dplyr::case_when(
        a_m2>1e9 & a_m2<1e10 ~ 2 # 100k ha to 1M ha
        , a_m2>=1e10 ~ max(round(a_m2*1e-11), 1)*3
        , T ~ 1
      )
      huge <- fact>1
    }else if (res>30){
      # 150B = 2; 200B = 2;...;500B = 5; 550B = 6;...; 900B = 9
      fact <- max(round(a_m2*1e-11), 1)
      huge <- fact>1
    }else{
      # this is a guess
      fact <- max(round(a_m2*1e-11), 1)*5
      huge <- fact>5
    }
    return(list(
      fact = fact
      , huge = huge
    ))
  }

  # # let's test this function
  # dplyr::tibble(
  #   a = seq(from = 0, to = 8e11, by = 0.5e9)
  # ) %>%
  # dplyr::rowwise() %>%
  # dplyr::mutate(
  #   fact = get_fact_fn(a, res=30) %>% purrr::pluck("fact")
  #   , a_ha = a/10000
  # ) %>%
  # dplyr::ungroup() %>%
  # ggplot(aes(x = a_ha, y = fact)) +
  #   geom_line() +
  #   labs(x = "area ha", y = "raster agg fact\nfor 30m rast") +
  #   scale_y_continuous(breaks = scales::extended_breaks(n=14)) +
  #   scale_x_log10(labels = scales::comma)

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! function 11
  # reclassify raster to set "undesirable" cells to NA so that these get filled with the nearest "desirable" value
  # foresttype "undesirable" = non-stocked, anything not in the lookup table
  # landfire "undesirable" = any CBD of <= 0, such as 0 and -999
  reclass_foresttype_rast <- function(rast, lookup){
    # mark all non-forest cells as "NA" in the cropped raster
      # rcl = two column matrix ("is", "becomes") can be useful for classifying integer
       # values. In that case, the arguments right and include.lowest are ignored.
      # unique type codes
      type_code_temp <- lookup %>%
        dplyr::pull(forest_type_code) %>%
        as.numeric() %>%
        unique()
    # matrix
      rcl_temp <- c(type_code_temp, type_code_temp) %>%
        matrix(ncol=2, byrow=F)
    # update raster to mark all non-forest cells as NA
      foresttype_temp <- rast %>%
        terra::classify(
          rcl = rcl_temp
          , others = NA
        )
    # return
      return(foresttype_temp)
  }
  reclass_landfire_rast <- function(rast){
    # only reclass if the raster is factor
    if(terra::is.factor(rast)){
      # mark all non-CBD cells as "NA" in the cropped raster
      # rcl = two column matrix ("is", "becomes") can be useful for classifying integer
       # values. In that case, the arguments right and include.lowest are ignored.
      lf_levels <- rast %>% terra::levels()
      rast <- as.numeric(rast)
      #divide by 100 to get actual kg/m^3
      lf_values_recode <- ifelse(
        lf_levels[[1]]$Value<=0
        , NA
        , lf_levels[[1]]$Value
        ) / 100
      # matrix
      rcl_temp <- c(lf_levels[[1]]$Value, lf_values_recode) %>%
        matrix(ncol=2, byrow=F)
      # update raster to mark all non-forest cells as NA and convert rest to numeric
      rast <- rast %>%
        terra::classify(
          rcl = rcl_temp
          , others = NA
        )
    }
    # return
    return(rast)
  }

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! function 14
  # we only need to use this function if we got points returned from the crop_raster_match_points()
  # where the values extracted from the raster were "NA" or "undesirable"
  # this function will check if we need to aggregate the raster data using get_fact_fn()
  # ...aggregate the raster data if needed using terra::aggregate()
  # ...fill in missing raster cells using fill_rast_na()
  ##### PARAMETERS
  # rast = raster returned from crop_raster_match_points()
  # ..... -OR- raster returned from crop_raster_match_points() which was then passed to:
  # ..... reclass_landfire_rast() or reclass_foresttype_rast()
  # bbox = bbox returned from crop_raster_match_points()
  # points = original spatial points passed to crop_raster_match_points()
  agg_fill_rast_match_points <- function(rast, bbox, points){
    ##################################
    # check if huge raster
    ##################################
      #  check if huge
      area_m2_temp <- bbox %>%
        sf::st_as_sfc() %>%
        terra::vect() %>%
        terra::project(terra::crs(rast)) %>% # this way we work with the same units no matter the input
        sf::st_as_sf() %>%
        sf::st_area() %>%
        as.numeric()

      # check it
      get_fact_fn_ans <- get_fact_fn(area_m2_temp, res = terra::res(rast)[1])
      is_huge_temp <- get_fact_fn_ans$huge

      # safe fill_rast_na first
      safe_fill_rast_na <- purrr::safely(fill_rast_na)

      ##### this process applies the fill_rast_na()
        # if is_huge_temp==T: aggregate, fill, update original rast resolution with filled
        # if is_huge_temp==F: fill
      if(is_huge_temp==T){
        # if huge, aggregate the cropped raster
        agg_rast <- rast %>%
          terra::aggregate(
            fact = get_fact_fn_ans$fact
            , fun = "modal", na.rm = T
            , cores = lasR::half_cores()
          )

        # `fill_rast_na()` process
        fill_rast_na_ans <- safe_fill_rast_na(agg_rast)

        # if no error keep going
        if(is.null(fill_rast_na_ans$error)){
            # get the filled result
            agg_rast <- fill_rast_na_ans$result
            # resample
            filler_rast_temp <- terra::resample(
              x = agg_rast
              , y = rast
              , method = "near"
            )
            # update the NA cells in the original raster
            rast <- terra::cover(
              x = rast
              , y = filler_rast_temp
              , values = NA
            )
        }

      }else{ # is_huge_temp==F
        # `fill_rast_na()` process
        fill_rast_na_ans <- safe_fill_rast_na(rast)

        # if no error keep going
        if(is.null(fill_rast_na_ans$error)){
            # get the filled result
            rast <- fill_rast_na_ans$result
        }

      }

    ##################################
    # if we successfully updated the raster
    # attach forest type to our original tree list
    ##################################
    if(terra::global(rast, fun = "isNA")==0){
      # now we apply the `terra::extract()` function to attach forest type to our original tree list
      tree_ftype_temp <-
        terra::extract(
          x = rast
          , y = points %>%
            terra::vect() %>%
            terra::project(terra::crs(rast)) # don't forget to reproject
        , ID = T # this way, the 2nd column will always be the raster value of the 1st layer
        ) %>%
        dplyr::rename(point_record_number = 1, raster_value = 2)
    }else{
      tree_ftype_temp <- NULL
    }

    # return
    return(list(
      point_values = tree_ftype_temp
      , rast = rast
    ))
  }
