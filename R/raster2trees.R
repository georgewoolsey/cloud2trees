#' @title Use a CHM raster to detect individual trees
#'
#' @description
#' `raster2trees()` is an all-in-one function to process a CHM raster
#' and return a spatial data frame of tree crown polygons and points.
#' The order of operations is:
#'
#' * Perform individual tree detection using [lidR::locate_trees()] with the [lidR::lmf()] algorithm
#' * Delineate tree crowns using [ForestTools::mcws()]
#'
#' Note, this function does not estimate DBH for the detected trees and only returns tree location, crown area, and height information.
#' To estimate tree DBH from the detected tree heights see [trees_dbh()].
#'
#' @param chm_rast raster. A  raster from `terra` or `stars`representing a canopy height model
#' @param outfolder string. The path of a folder to write the crown vector data to
#' @param ws numeric or function. Length or diameter of the moving window used to detect the local
#' maxima in the units of the input data (usually meters). If it is numeric a fixed window size is used.
#' If it is a function, the function determines the size of the window at any given location on the canopy.
#' By default function takes the height of a given pixel as its only argument and return the
#' desired size of the search window when centered on that pixel.
#' @param min_height numeric. Set the minimum height (m) for individual tree detection
#' @param min_crown_area numeric. Set the minimum crown area (m2) for individual tree detection
#' @param tempdir string. Directory to write intermediate files. Intermediate files are only created for large rasters too big to fit in memory.
#'
#' @references
#' [https://r-lidar.github.io/lidRbook/itd.html](https://r-lidar.github.io/lidRbook/itd.html)
#'
#' @return Returns a spatial data frame of individual tree crown vectors detected using the CHM.
#' The tree top point coordinates are located in the `tree_x` and `tree_y` columns.
#' The process also writes two `.gpkg` files to the `outfolder` directory: `chm_detected_crowns.gpkg` and `chm_detected_tree_tops.gpkg`
#'
#' @examples
#'  \dontrun{
#'  f <- paste0(system.file(package = "cloud2trees"),"/extdata/chm.tif")
#'  crowns_sf <- raster2trees(chm_rast = terra::rast(f), outfolder = tempdir())
#'  crowns_sf %>% class()
#'  crowns_sf %>% dplyr::glimpse()
#'  crowns_sf %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(fill=tree_height_m))
#'  }
#' @export
#'
raster2trees <- function(
  chm_rast
  , outfolder
  , ws = function(x) {
    y <- dplyr::case_when(
      is.na(x) ~ 1e-3 # requires non-null
      , x < 0 ~ 1e-3 # requires positive
      , x < 3.6 ~ 1.25 + x*0.15 # set lower bound
      , x > 32.5 ~ 5  # set upper bound
      , TRUE ~ exp( (0.0446*x) + (x^-0.555) ) # used gamma regression so exp the result
    )
    return(y)
  }
  , min_height = 2
  , min_crown_area = 0.1
  , tempdir = tempdir()
){
  ##############################################
  # check if raster is too big to fit in memory
  ##############################################
  split_raster <- split_raster_fn(chm_rast)

  ###############################################################################
  # make and process raster tiles if raster is too big to fit in memory
  ###############################################################################
  if(split_raster==T){
    # returns a list of files chm_tiles
    chm_tiles <- terra::makeTiles(
      chm_rast
      # specify rows and columns for tiles
      , y = c(
        round(terra::nrow(chm_rast)/4)
        , round(terra::ncol(chm_rast)/4)
      )
      , filename = file.path(tempdir,"tile_.tif")
      , na.rm = T
      , buffer = round(10/terra::res(chm_rast))[1] # 10m buffer
      , overwrite = T
    )
    # chm_tiles

    ##############################################
    # locate tree tops and crowns for tiles
    ##############################################
    # pass tiles to the function and return data frame of written files
    trees_crowns_data <- chm_tiles %>%
      purrr::map(\(x) trees_crowns_fn(
        rast_pth = x
        , ws = ws
        , hmin = min_height
        , dir = normalizePath(tempdir)
        , min_crown_area = min_crown_area
      )) %>%
      dplyr::bind_rows()

    write.csv(trees_crowns_data, file.path(tempdir, "trees_crowns_data.csv"), row.names = F)
    # trees_crowns_data = readr::read_csv(paste0(normalizePath(tempdir), "/trees_crowns_data.csv"))

    #################
    # bring tiles together
    #################
    # trees_crowns_data %>% head()
    # find overlapping buffers to select largest crown within overlaps
    rast_list_temp <- trees_crowns_data$chm_tile %>%
      purrr::map(function(x){
        r = terra::rast(x) %>%
          terra::classify(rcl = c(-Inf,Inf,1), others = 1)
        r[is.na(r)] = 1
        return(r)
      })

    # mosaic to get only buffers
    buffers_sf_temp <- terra::sprc(rast_list_temp) %>%
      terra::mosaic(fun = "sum") %>%
      terra::classify(rcl = c(2,Inf,1), others = NA) %>%
      # terra::plot(colNA="red")
      terra::as.polygons() %>%
      sf::st_as_sf() %>%
      sf::st_union()

    # read all crowns
    crowns_sf <- trees_crowns_data$crowns_sf_file %>%
      purrr::map(function(x){
        sf::st_read(x, quiet = T) %>%
          dplyr::mutate(fnm = basename(x))
      }) %>%
      dplyr::bind_rows()

    # get only crowns in buffer...make sure to return full crown and not just part in buffer
    buffer_crowns_temp <- crowns_sf %>%
      dplyr::inner_join(
        crowns_sf %>%
          sf::st_intersection(buffers_sf_temp) %>%
          sf::st_drop_geometry() %>%
          dplyr::distinct(layer,fnm)
        , by = dplyr::join_by(layer,fnm)
      )

    # remove buffer crowns from main file to add only selected crowns in later
    crowns_sf <- crowns_sf %>%
      dplyr::anti_join(
        buffer_crowns_temp %>% sf::st_drop_geometry()
        , by = dplyr::join_by(layer,fnm)
      )

    # list of tiles with crowns in buffered area
    fnm_list_temp <- buffer_crowns_temp$fnm %>% unique()

    # start with the first tile file
    if(nrow(buffer_crowns_temp)>0){
      keep_buffer_crowns_temp <- buffer_crowns_temp %>%
        dplyr::filter(fnm == fnm_list_temp[1])
    }else{keep_buffer_crowns_temp <- dplyr::slice_sample(crowns_sf, n = 0)}

    # sequentially add tiles and remove overlaps
    if(length(fnm_list_temp)>1){
      i <- 2
      while(i<=length(fnm_list_temp)){
        # read in next file
        x2_temp <- buffer_crowns_temp %>%
          dplyr::filter(fnm == fnm_list_temp[i])
        # remove overlaps
        keep_buffer_crowns_temp <- remove_overlap_fn(keep_buffer_crowns_temp, x2_temp, min_crown_area = min_crown_area)
        # clean and increment
        remove(x2_temp)
        i <- i+1
      }
    }

    # add keeps back to crowns data
    crowns_sf <- crowns_sf %>%
      dplyr::bind_rows(keep_buffer_crowns_temp) %>%
      # generate tree id
      dplyr::mutate(treeID = dplyr::row_number() %>% as.character()) %>%
      dplyr::relocate(treeID)

    # join tree tops
    tree_tops_temp <-
      1:nrow(trees_crowns_data) %>%
      purrr::map(function(x){
        sf::st_read(trees_crowns_data$tree_tops_file[x]) %>%
          dplyr::mutate(fnm = basename(trees_crowns_data$crowns_sf_file[x])) %>%
          dplyr::rename(
            layer = treeID
            , tree_height_m = Z
          ) %>%
          dplyr::mutate(
            tree_x = sf::st_coordinates(.)[,1]
            , tree_y = sf::st_coordinates(.)[,2]
          )
      }) %>%
      dplyr::bind_rows()
    # tree_tops_temp %>% dplyr::glimpse()

    # keep tree tops with crown match
      # note that the spatial location of a tree point might not fall within a crown vector
      # based on the vector trimming above...that's ok for our purposes as it won't impact the calculations
    crowns_sf <- crowns_sf %>%
      dplyr::inner_join(
        tree_tops_temp %>% sf::st_drop_geometry()
        , by = dplyr::join_by(layer==layer, fnm==fnm)
      ) %>%
      dplyr::select(-c(layer,fnm))
  }else{ # if(split_raster==T)
  ###############################################################################
  # if raster is not too big to fit in memory...it's simple
  ###############################################################################
    # identify tree tops
    tree_tops <- lidR::locate_trees(
      chm_rast
      , algorithm = lidR::lmf(
        ws = ws
        , hmin = min_height
      )
    )
    # check for validity
    if(
      !inherits(tree_tops, "sf") ||
      nrow(sf::st_zm(tree_tops, drop = T))==0 ||
      min(sf::st_is(tree_tops, type = c("POINT", "MULTIPOINT"))) == 0
    ){
      stop(paste0(
        "Could not locate any trees using the CHM raster"
        , "\n and window size settings...try different settings or data?"
        , "\n "
      ))
    }

    # delineate crowns
    crowns <- ForestTools::mcws(
      treetops = sf::st_zm(tree_tops, drop = T) # drops z values
      , CHM = chm_rast
      , minHeight = min_height
    )

    # ### Write the crown raster to the disk
    # terra::writeRaster(
    #   crowns
    #   , paste0(normalizePath(outfolder), "/top_down_detected_tree_crowns.tif")
    #   , overwrite = TRUE
    # )

    ### Convert crown raster to polygons, then to Sf
    crowns_sf <- crowns %>%
      # convert raster to polygons for each individual crown
      terra::as.polygons() %>%
      # fix polygon validity
      terra::makeValid() %>%
      # reduce the number of nodes in geometries
      terra::simplifyGeom() %>%
      # remove holes in polygons
      terra::fillHoles() %>%
      # convert to sf
      sf::st_as_sf() %>%
      dplyr::rename(layer = 1) %>%
      # get the crown area
      dplyr::mutate(
        crown_area_m2 = as.numeric(sf::st_area(.))
      ) %>%
      #remove super small crowns
      dplyr::filter(
        crown_area_m2 > min_crown_area
      )

      ### add tree data to tree_tops
      tree_tops = tree_tops %>%
        # pull out the coordinates and update treeID
        dplyr::mutate(
          tree_x = sf::st_coordinates(.)[,1]
          , tree_y = sf::st_coordinates(.)[,2]
          , tree_height_m = sf::st_coordinates(.)[,3]
          , treeID = paste(treeID,round(tree_x, 1),round(tree_y, 1), sep = "_")
        )
      # str(tree_tops)

      ### Join the crowns with the tree tops to append data, remove Nulls
      crowns_sf = crowns_sf %>%
        sf::st_join(tree_tops) %>%
        dplyr::group_by(layer) %>%
        dplyr::mutate(n_trees = dplyr::n()) %>%
        dplyr::group_by(treeID) %>%
        dplyr::mutate(n_crowns = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::filter(
          n_trees==1
          | (n_trees>1 & n_crowns==1) # keeps the treeID that only has one crown if multiple trees to one crown
        ) %>%
        dplyr::select(c(
          "treeID", "tree_height_m"
          , "tree_x", "tree_y"
          , "crown_area_m2"
          , tidyselect::starts_with("comp_")
        )) %>%
        dplyr::filter(!is.na(treeID))
      # str(crowns_sf)

      # ## !!!!!!!!!!!!!!!!!!!!!!! Not needed...and doesn't work with tiled CHM process
      # ### Add crown data summaries
      # crown_sum_temp = data.frame(
      #     mean_crown_ht_m = terra::extract(x = chm_rast, y = terra::vect(crowns_sf), fun = "mean", na.rm = T)
      #     , median_crown_ht_m = terra::extract(x = chm_rast, y = terra::vect(crowns_sf), fun = "median", na.rm = T)
      #     , min_crown_ht_m = terra::extract(x = chm_rast, y = terra::vect(crowns_sf), fun = "min", na.rm = T)
      #   ) %>%
      #   dplyr::select(-c(tidyselect::ends_with(".ID"))) %>%
      #   dplyr::rename_with(~ stringr::str_remove_all(.x,".Z")) %>%
      #   dplyr::rename_with(~ stringr::str_remove_all(.x,".focal_mean"))
      # ### join crown data summary
      # crowns_sf = crowns_sf %>%
      #   dplyr::bind_cols(crown_sum_temp)
      # # str(crowns_sf)
  }
  ##############################################
  ### write the data to the disk
  ##############################################
  if(nrow(crowns_sf)>250e3){
      # split up the detected crowns
      crowns_sf = crowns_sf %>%
        dplyr::arrange(as.numeric(tree_x),as.numeric(tree_y)) %>%
        # groups of 500k
        dplyr::mutate(grp = ceiling(dplyr::row_number()/250e3))

      write_fnl_temp = crowns_sf$grp %>%
        unique() %>%
        purrr::map(function(x){
          ### write the data to the disk
          # crown vector polygons
          sf::st_write(
            crowns_sf %>%
              dplyr::filter(grp == x) %>%
              dplyr::select(-c(grp))
            , paste0(normalizePath(outfolder), "/chm_detected_crowns_",x,".gpkg")
            , append = FALSE
            , quiet = TRUE
          )
          # tree top vector points
          sf::st_write(
            # get tree points
            crowns_sf %>%
              dplyr::filter(grp == x) %>%
              dplyr::select(-c(grp)) %>%
              sf::st_drop_geometry() %>%
              sf::st_as_sf(coords = c("tree_x", "tree_y"), crs = sf::st_crs(crowns_sf))
            , paste0(normalizePath(outfolder), "/chm_detected_tree_tops_",x,".gpkg")
            , append = FALSE
            , quiet = TRUE
          )
          return(
            dplyr::tibble(
              crowns_file = paste0(normalizePath(outfolder), "/chm_detected_crowns_",x,".gpkg")
              , trees_file = paste0(normalizePath(outfolder), "/chm_detected_tree_tops_",x,".gpkg")
            )
          )
        }) %>%
        dplyr::bind_rows()
    }else{
        # crown vector polygons
        sf::st_write(
          crowns_sf
          , paste0(normalizePath(outfolder), "/chm_detected_crowns.gpkg")
          , append = FALSE
          , quiet = TRUE
        )
        # tree top vector points
        sf::st_write(
          # get tree points
          crowns_sf %>%
            sf::st_drop_geometry() %>%
            sf::st_as_sf(coords = c("tree_x", "tree_y"), crs = sf::st_crs(crowns_sf))
          , paste0(normalizePath(outfolder), "/chm_detected_tree_tops.gpkg")
          , append = FALSE
          , quiet = TRUE
        )
    }
  # return
  return(crowns_sf)
}

###___________________________________________________###
### Intermediate functions
###___________________________________________________###

###____________________________________###
### 1/3 compare 2 crown tiles and
### get rid of overlaps by selecting
### the largest crown for each overlap
###____________________________________###
  remove_overlap_fn <- function(x1, x2, min_crown_area){
    # make features valid
    # ... without this: Error: TopologyException: Input geom 0 is invalid: Nested shells
    # ... https://github.com/r-spatial/sf/issues/870
      # primary data
      x1 <- x1 %>%
        sf::st_make_valid() %>%
        dplyr::filter(
          sf::st_is_valid(.)
          & !sf::st_is_empty(.)
        )
      if(nrow(x1)==0){return(x2)} ## could issue error here instead
      # data to check overlap
      x2 <- x2 %>%
        sf::st_make_valid() %>%
        dplyr::filter(
          sf::st_is_valid(.)
          & !sf::st_is_empty(.)
        )
      if(nrow(x2)==0){return(x1)}

    # identify equal vectors
      equals_temp <- x1 %>%
        sf::st_join(x2, join = sf::st_equals, left = F) %>%
        sf::st_drop_geometry()

    # get spatial intersection and filter
      intersect_temp <- x1 %>%
        sf::st_intersection(x2) %>%
        dplyr::mutate(area = sf::st_area(.) %>% as.numeric) %>%
        sf::st_drop_geometry() %>%
        dplyr::filter(area>0) %>%
        # remove exact matches in x1
        dplyr::anti_join(
          equals_temp
          , by = dplyr::join_by(layer == layer.x, fnm == fnm.x)
        ) %>%
        # remove exact matches in x2
        dplyr::anti_join(
          equals_temp
          , by = dplyr::join_by(layer.1 == layer.y, fnm.1 == fnm.y)
        ) %>%
        # make one comparison based on x2...similar to sf::st_join(largest = T)
        dplyr::group_by(layer.1,fnm.1) %>%
        dplyr::filter(area == max(area)) %>%
        dplyr::ungroup() %>%
        # what NOT to keep...get list of polygons that intersect and get rid of the smaller one
        dplyr::mutate(
          layer = dplyr::case_when(
            crown_area_m2 > crown_area_m2.1 ~ layer.1
            , crown_area_m2 < crown_area_m2.1 ~ layer
            , T ~ layer.1
          )
          , fnm = dplyr::case_when(
            crown_area_m2 > crown_area_m2.1 ~ fnm.1
            , crown_area_m2 < crown_area_m2.1 ~ fnm
            , T ~ fnm.1
          )
        ) %>%
        dplyr::distinct(layer, fnm)

    x1_new <- x1 %>%
      # filter to keep biggest when area is not equal
      dplyr::anti_join(
        intersect_temp
        , by = dplyr::join_by("layer","fnm")
      )

    x2_new <- x2 %>%
      # filter to keep biggest when area is not equal
      dplyr::anti_join(
        intersect_temp
        , by = dplyr::join_by("layer","fnm")
      ) %>%
      # filter to keep x1 when area is equal
      dplyr::anti_join(
        equals_temp %>%
          dplyr::distinct(layer.y, fnm.y)
        , by = dplyr::join_by(layer == layer.y, fnm == fnm.y)
      )

    # identify remaining overlaps
    olap_temp <- x1_new %>%
      sf::st_intersection(x2_new) %>%
      dplyr::mutate(area = sf::st_area(.) %>% as.numeric()) %>%
      dplyr::filter(area > 0) %>%
      sf::st_union()

    if(dplyr::coalesce(length(olap_temp),0)>0){
      # combine filtered crowns
      new_dta <- dplyr::bind_rows(x1_new, x2_new) %>%
        # remove remaining overlaps
        sf::st_difference(olap_temp) %>%
        dplyr::mutate(area = sf::st_area(.) %>% as.numeric()) %>%
        sf::st_make_valid() %>%
        dplyr::filter(
          area > min_crown_area
          & sf::st_is_valid(.)
          & !sf::st_is_empty(.)
        ) %>%
        dplyr::select(-c(area))
    }else{
      # combine filtered crowns
      new_dta <- dplyr::bind_rows(x1_new, x2_new) %>%
        dplyr::mutate(area = sf::st_area(.) %>% as.numeric()) %>%
        sf::st_make_valid() %>%
        dplyr::filter(
          area > min_crown_area
          & sf::st_is_valid(.)
          & !sf::st_is_empty(.)
        ) %>%
        dplyr::select(-c(area))
    }

    return(new_dta)
  }

###____________________________________###
### 2/3 locate tree tops and crowns
### for tiled rasters
###____________________________________###
  #################
  # function to locate tree tops and crowns for each tile
  #################
  trees_crowns_fn <- function(rast_pth, ws, hmin, dir, min_crown_area){
    # read raster
    r <- terra::rast(rast_pth)
    ###################
    # locate seeds
    ###################
    tree_tops <- lidR::locate_trees(
      r
      , algorithm = lidR::lmf(
        ws = ws
        , hmin = hmin
      )
    )

    # write
    tree_tops_file <- paste0(
      dir
      ,"/tree_tops_"
      , stringr::str_replace_all(basename(rast_pth), pattern = ".tif", replacement = ".gpkg")
    )

    sf::st_write(
      tree_tops
      , tree_tops_file
      , append = F
      , quiet = T
    )
    ###################
    # crowns
    ###################
    crowns = ForestTools::mcws(
      treetops = sf::st_zm(tree_tops, drop = T) # drops z values
      , CHM = r
      , minHeight = hmin
    )

    # write
    crowns_file = paste0(
      dir
      ,"/crowns_"
      , basename(rast_pth)
    )

    terra::writeRaster(
      crowns
      , filename = crowns_file
      , overwrite = T

    )

    ###################
    # crowns vector
    ###################
    crowns_sf = crowns %>%
      # convert raster to polygons for each individual crown
      terra::as.polygons() %>%
      # fix polygon validity
      terra::makeValid() %>%
      # reduce the number of nodes in geometries
      terra::simplifyGeom() %>%
      # remove holes in polygons
      terra::fillHoles() %>%
      # convert to sf
      sf::st_as_sf() %>%
      dplyr::rename(layer = 1) %>%
      # get the crown area
      dplyr::mutate(
        crown_area_m2 = as.numeric(sf::st_area(.))
      ) %>%
      #remove super small crowns
      dplyr::filter(
        crown_area_m2 > min_crown_area
      )

    # write
    crowns_sf_file = stringr::str_replace_all(crowns_file, pattern = ".tif", replacement = ".gpkg")

    sf::st_write(
      crowns_sf
      , crowns_sf_file
      , append = F
    )

    # return
    return(
      dplyr::tibble(
        chm_tile = rast_pth
        , tree_tops_file = tree_tops_file
        , crowns_terra_file = crowns_file
        , crowns_sf_file = crowns_sf_file
      )
    )
  }
###____________________________________###
### 3/3 check if split raster for
### large on-disk rasters
###____________________________________###
  # # !!!!!! lidR::locate_trees didn't work with large rasters, given error:
  # # !!!!!! Error: Large on-disk rasters are not supported by locate_tree. Load the raster manually.
  # # ... see: https://github.com/r-lidar/lidR/blob/9e052574adb40513e6d88ac141cfd0d41c3798ba/R/utils_raster.R
  # # ... see: https://stackoverflow.com/questions/77863431/manually-load-large-raster-to-perform-lidrlocate-trees
  # # check if need to split raster into tiles
  split_raster_fn <- function(raster, n = 10){
    # puts raster in memory = raster_in_memory
    raster <- raster*1
    # check if is not in memory = raster_is_proxy
    raster_is_proxy <- !terra::inMemory(raster)
    # check if memory is available = raster_fits_in_memory
    if(raster_is_proxy){
      # check memory needed
      nc <- terra::nrow(raster)*terra::ncol(raster)
      n <- n*terra::nlyr(raster)
      memneed <- nc * n * 8L
      memavail <- terra::free_RAM()*1000
      memavail <- 0.6 * memavail
      is_mem_avail <- memneed < memavail
    }else{is_mem_avail <- T}
    # return check for splitting up raster into tiles
    split <- dplyr::case_when(
      raster_is_proxy == F ~ F
      , raster_is_proxy == T & is_mem_avail == T ~ F
      , raster_is_proxy == T & is_mem_avail == F ~ T
    )
    return(split)
  }
