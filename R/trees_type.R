#' @title Estimate forest type for a tree list based on location
#'
#' @description
#' `trees_type()` uses the input tree list (e.g. as exported by [raster2trees()]) with the columns
#' `treeID`, `tree_x`, `tree_y` to attach species information using USDA Forest Inventory and Analysis (FIA) codes.
#' If a spatial data frame of points is the input tree list, then the columns `tree_x`, `tree_y` are not required.
#'
#' FIA Forest Type Group Code is attached to each tree in the tree list based on the spatial overlap with the Forest Type Groups of the Continental United States dataset [Wilson 2023](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d).
#'
#' The simplified process for attaching forest type group to a tree is:
#'
#' * Forest type group 30-m raster (Wilson 2023) was aggregated to 90-m to make the data more accessible over the entire continental US
#' * Nearest neighbor imputation is used to fill forest type data if a tree falls inside a no-forest cell in the original data
#' * The FIA forest type group is applied to a tree based on spatial overlap
#'
#' @param tree_list data.frame. A data frame with the columns `treeID`, `tree_x`, `tree_y`, and `tree_height_m`.
#' If an `sf` class object with POINT geometry (see [sf::st_geometry_type()]), the program will use the data "as-is" and only require the `treeID` column.
#' @param crs string. A crs string as returned from [sf::st_crs()] or the EPSG code of the x,y coordinates.
#' Defaults to the crs of the `tree_list` data if of class "sf".
#' @param study_boundary sf. The boundary of the study are to define the area of the regional model.
#' If no boundary given, regional model will be built from location of trees in the tree list.
#' @param input_foresttype_dir directory where Forest Type Groups data exists. Use [get_foresttype()] first.
#' @param max_search_dist_m number. Maximum search distance (m) to obtain forest type group data for trees in `tree_list` that overlap with non-forest data in the original Wilson (2023) data.
#' Larger search distances will increase processing time and possibly result in memory issues.
#'
#' @references
#' * [Forest Type Groups of the Continental United States](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d)
#' Wilson, B.T. (2023). Forest Type Groups of the Continental United States.
#'
#' @return Returns a list of objects: tree_list = spatial data frame of individual trees; foresttype_rast = raster of forest types in the area.
#'
#' @examples
#'  \dontrun{
#'  # example tree list
#'  tl <- dplyr::tibble(
#'      treeID = c(1:21)
#'      , tree_x = rnorm(n=21, mean = 458064, sd = 11)
#'      , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
#'    )
#'  # call the function
#'  tl_type <- trees_type(tree_list = tl, crs = "32613")
#'  # what?
#'  tl_type %>% class()
#'  # a list, but what is in it?
#'  tl_type %>% names()
#'  # plot the tree_list spatial points
#'  tl_type$tree_list %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=forest_type_group))
#'  # plot the foresttype_rast raster
#'  tl_type$foresttype_rast %>% terra::plot()
#'  }
#' @export
#'
trees_type <- function(
  tree_list
  , crs = NA
  , study_boundary = NA
  , input_foresttype_dir = NULL
  , max_search_dist_m = 1000
){
  ####################################################################
  # check external data
  ####################################################################
    # find external data
    find_ext_data_ans <- find_ext_data(
      input_foresttype_dir = input_foresttype_dir
    )
    # if can't find external foresttype data
    if(is.null(find_ext_data_ans$foresttype_dir)){
      stop(paste0(
        "Forest Type Group data has not been downloaded to package contents. Use `get_foresttype()` first."
        , "\nIf you supplied a value to the `input_foresttype_dir` parameter check that directory for data."
      ))
    }
  ##################################
  # ensure that tree data exists
  ##################################
  f <- tree_list %>% names()
  if(length(f)==0){f <- ""}
  if(
    max(grepl("treeID", f))==0
  ){
    stop(paste0(
      "`tree_list` data must contain `treeID` column."
      , "\nProvide the `treeID` as a unique identifier of individual trees."
    ))
  }

  ##################################
  # convert to spatial points data
  ##################################
  if(inherits(tree_list, "sf")){
    # if points, just use it
    if( min(sf::st_is(tree_list, type = c("POINT", "MULTIPOINT"))) == 1 ){
      tree_tops <- tree_list %>%
        dplyr::mutate(
          tree_x = sf::st_coordinates(.)[,1]
          , tree_y = sf::st_coordinates(.)[,2]
        )
    }else{ # if spatial but not points, drop geom and set to points
      if(
        max(grepl("tree_x", names(tree_list)))==0
        || max(grepl("tree_y", names(tree_list)))==0
      ){ # doesn't contain x,y
        tree_tops <- tree_list %>%
          sf::st_centroid() %>%
          dplyr::mutate(
            tree_x = sf::st_coordinates(.)[,1]
            , tree_y = sf::st_coordinates(.)[,2]
          )
      }
      tree_tops <- tree_list %>%
        sf::st_drop_geometry() %>%
        sf::st_as_sf(
          coords = c("tree_x", "tree_y"), crs = sf::st_crs(tree_list)
          , remove = F
        ) %>%
        dplyr::mutate(treeID = as.character(treeID))
    }
  }else{ # not spatial data
    # convert from data.frame to spatial points
    if(!inherits(tree_list, "data.frame")){
      stop("must pass a data.frame or sf object to the `tree_list` parameter")
    }
    if(is.na(crs) || is.na(readr::parse_number(as.character(crs)))){
      stop("must provide the EPSG code in `crs` parameter for the projection of x,y data")
    }
    tree_tops <- tree_list %>%
      sf::st_as_sf(
        coords = c("tree_x", "tree_y")
        , crs = paste0( "EPSG:", readr::parse_number(as.character(crs)) )
        , remove = F
      ) %>%
      dplyr::mutate(treeID = as.character(treeID))
  }

  # check for duplicate treeID
  if(
    nrow(tree_tops) != length(unique(tree_tops$treeID))
  ){
    stop("Duplicates found in the treeID column. Please remove duplicates and try again.")
  }
  # get rid of columns we'll create
    tree_tops <- tree_tops %>%
      # throw in hey_xxxxxxxxxx to test it works if we include non-existant columns
      dplyr::select( -dplyr::any_of(c(
        "hey_xxxxxxxxxx"
        , "forest_type_group_code"
        , "forest_type_group"
        , "hardwood_softwood"
      )))

  ##################################
  # load foresttype data. see get_foresttype()
  ##################################
    # get the foresttype data
    foresttype <- terra::rast(
        file.path(find_ext_data_ans$foresttype_dir, "foresttype.tif")
      )
    # get the lookup
    foresttype_lookup <- readr::read_csv(
        file.path(find_ext_data_ans$foresttype_dir, "foresttype_lookup.csv")
        , progress = F
        , show_col_types = F
      ) %>%
      dplyr::mutate(dplyr::across(
        dplyr::everything()
        , as.character
      ))

  ##################################
  # define extent to crop forest type raster
  ##################################
    # get extent of trees data
    bbox_temp <- tree_tops %>% sf::st_bbox()
    # find largest side
    buffer_temp <- max(
        bbox_temp["xmax"]-bbox_temp["xmin"]
        , bbox_temp["ymax"]-bbox_temp["ymin"]
      ) %>%
      max(max_search_dist_m)

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
        )
      # reset bbox and buffer if larger
      if(buffer_b_temp>buffer_temp){
        buffer_temp <- buffer_b_temp
        bbox_temp <- bbox_b_temp
      }
    }

    # apply the buffer to get the extent
    ext_temp <- bbox_temp %>%
      sf::st_as_sfc() %>%
      sf::st_buffer(buffer_temp, endCapStyle = "SQUARE") %>%
      terra::vect() %>%
      terra::project(terra::crs(foresttype))

  ##################################
  # mask the raster to this extent
  ##################################
    foresttype_temp <- foresttype %>%
      terra::crop(ext_temp) # crop it first to make if faster

  ##################################
  # do we even need to fill NA values for this list?
  ##################################
    # attempt to extract the forest type for the trees
    # now we apply the `terra::extract()` function to attach forest type to our original tree list
    tree_ftype_temp <-
      terra::extract(
        x = foresttype_temp
        , y = tree_tops %>%
          terra::vect() %>%
          terra::project(terra::crs(foresttype)) # don't forget to reproject
      )
    # let's check with the lookup table
    na_trees <- tree_ftype_temp %>%
      dplyr::mutate(
        forest_type_code = foresttype %>% as.character()
      ) %>%
      dplyr::left_join(foresttype_lookup, by = "forest_type_code") %>%
      dplyr::filter(is.na(forest_type_group_code)) %>%
      nrow() %>%
      dplyr::coalesce(1)
  ##################################
  # IF we even need to fill NA values
  ##################################
  if(na_trees>0){
    # mark all non-forest cells as "NA" in the cropped raster
      # rcl = two column matrix ("is", "becomes") can be useful for classifying integer
       # values. In that case, the arguments right and include.lowest are ignored.
      # unique type codes
      type_code_temp <- foresttype_lookup %>%
        dplyr::pull(forest_type_code) %>%
        as.numeric() %>%
        unique()
    # matrix
      rcl_temp <- c(type_code_temp, type_code_temp) %>%
        matrix(ncol=2, byrow=F)
    # update raster to mark all non-forest cells as NA
      foresttype_temp <- foresttype_temp %>%
        terra::classify(
          rcl = rcl_temp
          , others = NA
        )

    ##################################
    # check if huge raster
    ##################################
      #  check if huge
      area_m2_temp <- bbox_temp %>%
        sf::st_as_sfc() %>%
        terra::vect() %>%
        terra::project(terra::crs(foresttype)) %>% # this way we work with the same units no matter the input
        sf::st_as_sf() %>%
        sf::st_area() %>%
        as.numeric()

      # check it
      get_fact_fn_ans <- get_fact_fn(area_m2_temp, res = terra::res(foresttype_temp)[1])
      is_huge_temp <- get_fact_fn_ans$huge

      # safe fill_rast_na first
      safe_fill_rast_na <- purrr::safely(fill_rast_na)

      ##### this process applies the fill_rast_na()
        # if is_huge_temp==T: aggregate, fill, update original rast resolution with filled
        # if is_huge_temp==F: fill
      if(is_huge_temp==T){
        # if huge, aggregate the cropped raster
        agg_foresttype_temp <- foresttype_temp %>%
          terra::aggregate(
            fact = get_fact_fn_ans$fact
            , fun = "modal", na.rm = T
            , cores = lasR::half_cores()
          )

        # `fill_rast_na()` process
        fill_rast_na_ans <- safe_fill_rast_na(agg_foresttype_temp)

        # if no error keep going
        if(is.null(fill_rast_na_ans$error)){
            # get the filled result
            agg_foresttype_temp <- fill_rast_na_ans$result
            # resample
            filler_rast_temp <- terra::resample(
              x = agg_foresttype_temp
              , y = foresttype_temp
              , method = "near"
            )
            # update the NA cells in the original raster
            foresttype_temp <- terra::cover(
              x = foresttype_temp
              , y = filler_rast_temp
              , values = NA
            )
        }

      }else{
        # `fill_rast_na()` process
        fill_rast_na_ans <- safe_fill_rast_na(foresttype_temp)

        # if no error keep going
        if(is.null(fill_rast_na_ans$error)){
            # get the filled result
            foresttype_temp <- fill_rast_na_ans$result
        }

      }

    ##################################
    # if we successfully updated the raster
    # attach forest type to our original tree list
    ##################################
    if(terra::global(foresttype_temp, fun = "isNA")==0){
      # now we apply the `terra::extract()` function to attach forest type to our original tree list
      tree_ftype_temp <-
        terra::extract(
          x = foresttype_temp
          , y = tree_tops %>%
            terra::vect() %>%
            terra::project(terra::crs(foresttype)) # don't forget to reproject
        )
    }

  }
  #######################################################
  # prep final data
  #######################################################
    if(length(tree_ftype_temp$foresttype)==nrow(tree_tops)){
      # now let's join it back with our data and check it
      tree_tops <- tree_tops %>%
        dplyr::mutate(
          forest_type_code = tree_ftype_temp$foresttype %>% as.character()
        ) %>%
        dplyr::left_join(
          foresttype_lookup %>%
            dplyr::select(
              forest_type_code
              , forest_type_group_code
              , forest_type_group
              , hardwood_softwood
            )
          , by = "forest_type_code"
        ) %>%
        dplyr::select(-forest_type_code)
    }else{ # if
      tree_tops <- tree_tops %>%
        dplyr::mutate(
          forest_type_group_code = as.character(NA)
          , forest_type_group = as.character(NA)
          , hardwood_softwood = as.character(NA)
        )
      message(paste0(
        "Unable to determine forest type for this tree list and study boundary (if provided)."
        , "\nTry expanding the study boundary area or increasing the max_search_dist_m parameter"
        , "\nand ensure that your tree data is in the continental US."
      ))
    }

  # return
  return(list(
    tree_list = tree_tops
    , foresttype_rast = foresttype_temp
  ))
}
#######################################################
# intermediate function 1
#######################################################
  # fill NA values of raster using nearest neighbor interpolation
  fill_rast_na <- function(rast){
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
      fact <- max(round(a_m2*1e-11), 1)*3
      huge <- fact>3
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
