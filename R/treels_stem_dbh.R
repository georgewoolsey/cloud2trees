#' @title detect tree stems and estimate DBH using `TreeLS` package
#'
#' @description
#' `treels_stem_dbh()` is an all-in-one function to process height normalized .las|.laz files
#' using the functionality of the `TreeLS` package. The function generates a list of stems and estimates DBH directly from the point cloud.
#' `treels_stem_dbh()` outputs:
#'
#' * A .laz file with the `Classification` data updated to: ground points (class 2); water points (class 9); stem points (class 4); non-stem (class 5).
#' * A vector data file in `gpkg` format with the tree identification stem locations, heights, and DBH estimates.
#'
#' The order of operations is:
#'
#' * Detect tree stems/boles from the height normalized point cloud using [TreeLS::treeMap()] with the [TreeLS::map.hough()]  algorithm
#' * Merge overlapping tree coordinates using [TreeLS::treeMap.merge()]
#' * Assign tree IDs to the original points using [TreeLS::treePoints()] with the [TreeLS::trp.crop()] algorithm
#' * Flag only the stem points using [TreeLS::stemPoints()] with the [TreeLS::stm.hough()] algorithm
#' * Perform DBH estimation using [TreeLS::tlsInventory()] with the [TreeLS::shapeFit()] algorithm
#'
#' @param folder string. The path of a folder containing a set of las/laz files. Can also be a vector of file paths.
#' @param outfolder string. The path of a folder to write the tiled las and vector files to
#' @param min_height numeric. Set the minimum height (m) for individual tree detection
#' @param max_dbh numeric. Set the largest tree diameter (m) expected in the point cloud
#' @param chunk_these logical. Do the las/laz files need to be tiled to work with smaller subsets? See `is_chunked_grid` in [chunk_las_catalog()]
#'
#' @references
#' https://github.com/tiagodc/TreeLS
#'
#' @return Returns an `sf` data.frame with TreeLS detected trees and DBH estimated directly from the point cloud.
#' Exports files in the
#' `outfolder` defined by the user in the function call.
#'
#' @examples
#'  \dontrun{
#'  o <- "../data"
#'  i <- "../data/normlasdata"
#'  r <- cloud2trees::treels_stem_dbh(folder = i, outfolder = o)
#'  r %>% names()
#'  }
#' @export
#'
treels_stem_dbh <- function(
    folder
    , outfolder
    , min_height = 2
    , max_dbh = 2
    , chunk_these = FALSE
) {
  # chunk las files with buffer if needed
  if(chunk_these==T){
    normalize_flist <- chunk_norm_las(folder = folder, outfolder = outfolder)
  }else{ # otherwise get list of files to process
    # get list of files to process
    # check if folder contains las files directly
    chk <- folder %>%
      tolower() %>%
      stringr::str_detect(".*\\.(laz|las)$")
    # file list
    if(max(chk)==1){
      ff <- stringr::str_subset(folder, ".*\\.(laz|las)$")
    }else{ # otherwise just search the folder
      ff <- list.files(normalizePath(folder), pattern = ".*\\.(laz|las)$", full.names = T)
    }
    # create spatial index files (.lax)
    normalize_flist <- create_lax_for_tiles(las_file_list = ff)
  }

  # map over the normalized point cloud tiles
  treels_pipeline_ans <-
    normalize_flist %>%
    purrr::map(\(x) treels_pipeline(
      las = x
      , outfolder = outfolder
      , min_height = min_height
      , max_dbh = max_dbh
    ))

  # create spatial index files (.lax)
  las_stem_flist <- create_lax_for_tiles(
    list.files(normalizePath(outfolder), pattern = ".*\\.(laz|las)$", full.names = T)
  )

  # Combine the vector data written as `gpkg` tile files
  ###__________________________________________________________###
  ### Merge the stem vector location tiles into a single object ###
  ###__________________________________________________________###
  if(
    length(list.files(normalizePath(outfolder), pattern = ".*\\.gpkg$", full.names = T)) > 0
  ){
    treels_dbh_locations <- list.files(normalizePath(outfolder), pattern = ".*\\.gpkg$", full.names = T) %>%
        purrr::map(\(x) sf::st_read(
          dsn = x
          , quiet = T
        )) %>%
        dplyr::bind_rows() %>%
        sf::st_as_sf() %>%
        sf::st_make_valid()
        # sf::st_set_crs(proj_crs)

    ## Clean the Stem Vector Data
    # The cleaning process uses the following steps:
    # * remove stems with empty radius estimates from the `TreeLS::tlsInventory` DBH estimation step
    # * remove stems >= DBH threshold set by the user in the parameter `max_dbh` (`r max_dbh`m in this example)
    # * remove stems with empty or invalid xy coordinates

    treels_dbh_locations <- treels_dbh_locations %>%
      dplyr::filter(
        !is.na(radius_m)
        & dbh_m <= max_dbh
        & sf::st_is_valid(.)
        & !sf::st_is_empty(.)
      ) %>%
      dplyr::mutate(
        condition = "detected_stem"
      )

    # ###___________________________________________________________###
    # ### Write the detected DBHs
    # ###___________________________________________________________###
    #   sf:::st_write(
    #     treels_dbh_locations
    #     , dsn = paste0(config$delivery_dir, "/bottom_up_detected_stem_locations.gpkg")
    #     , append = FALSE
    #     , delete_dsn = TRUE
    #     , quiet = TRUE
    #   )
  }else{treels_dbh_locations <- dplyr::tibble(NULL)}

  # return
  return(
    treels_dbh_locations
  )
}
###_____________________________________________________###
### Intermediate functions ###
###_____________________________________________________###

###___________________________________________###
### 1/3 Define function to map for potential  ###
### tree locations w TreeLS::treeMap          ###
###___________________________________________###
  ### Function to map for potential tree locations with error handling
  treels_treemap <- function(las, max_dbh = 2){
    result <- tryCatch(
      expr = {
        map = TreeLS::treeMap(
          las = las
          , method = TreeLS::map.hough(
            # height thresholds applied to filter a point cloud before processing
            # this is for detecting stems...not determining tree height
            min_h = 1
            , max_h = 5
            # height interval to perform point filtering/assignment/classification
            , h_step = 0.5
            # pixel side length to discretize the point cloud layers
              # while performing the Hough Transform circle search
            , pixel_size = 0.025
            # largest tree diameter expected in the point cloud
            , max_d = max_dbh # 0.75m = 30in
            # minimum point density (0 to 1) within a pixel evaluated
              # on the Hough Transform - i.e. only dense point clousters will undergo circle search
              # hey google, define "clouster" ?
            , min_density = 0.0001
            # minimum number of circle intersections over a pixel
              # to assign it as a circle center candidate.
            , min_votes = 3
          )
          # parameter passed down to treeMap.merge (if merge > 0)
          , merge = 0
        )
      },
      error = function(e) {
        message <- paste("Error:", e$message)
        return(message)
      }
    )
    if (inherits(result, "error")) {
      return(result)
    } else {
      return(result)
    }
  }
###___________________________________________###
### 2/3 input single norm las then use        ###
### TreeLS flow to classify+dbh               ###
###___________________________________________###
  # pass this function a file path of the normalized las you wish to detect stems and classify
  treels_pipeline <- function(las, outfolder, min_height = 2, max_dbh = 2) {
      ### Get the desired las file
      las_name <- basename(las)

      ### Read in the desired las file
      las_norm_tile <- lidR::readLAS(las)
      las_norm_tile <- lidR::filter_poi(las_norm_tile, Z >= 0)

      # get the maximum point height
      max_point_height <- max(las_norm_tile@data$Z)

      # IF MAX HEIGHT GOOD...KEEP DOING IT
      if(max_point_height >= min_height){
        ###______________________________________________________________###
        ### 1) Apply the `TreeLS::treeMap` [stem detection function](#detect_stem_fn)
        ###______________________________________________________________###
        ### Run the function to search for candidate locations
        treemap_temp <- treels_treemap(las_norm_tile, max_dbh = max_dbh)

        ### If the class of the result == "LAS"...REALLY KEEP DOING IT
        if(inherits(treemap_temp, "LAS")){
          ###______________________________________________________________###
          ### 2) Merge overlapping tree coordinates using `TreeLS::treeMap.merge`
          ###______________________________________________________________###
          treemap_temp <- TreeLS::treeMap.merge(treemap_temp)
          ###______________________________________________________________###
          ### 3) Assign tree IDs to the original points using `TreeLS::treePoints`
          ###______________________________________________________________###
          ### Classify tree regions
          ## Assigns TreeIDs to a LAS object based on coordinates extracted from a treeMap object.
          las_norm_tile <- TreeLS::treePoints(
            las = las_norm_tile
            , map = treemap_temp
            , method = TreeLS::trp.crop(l = 3)
          )
          # plot(las_norm_tile, color = "TreeID")

          ###______________________________________________________________###
          ### 4) Flag only the stem points using `TreeLS::stemPoints`
          ###______________________________________________________________###
          ### Classify stem points
          las_norm_tile <- TreeLS::stemPoints(
            las = las_norm_tile
            , method = TreeLS::stm.hough(
              # height interval to perform point filtering/assignment/classification.
              h_step = 0.5
              # largest tree diameter expected in the point cloud
              , max_d = max_dbh # 0.75m = 30in
              # tree base height interval to initiate circle search
              , h_base = c(1, 2.5)
              #  pixel side length to discretize the point cloud layers
                # while performing the Hough Transform circle search.
              , pixel_size = 0.025
              # minimum point density (0 to 1) within a pixel evaluated
                # on the Hough Transform - i.e. only dense point clousters will undergo circle search
                # hey google, define "clouster" ?
              , min_density = 0.1
              # minimum number of circle intersections over a pixel
                # to assign it as a circle center candidate.
              , min_votes = 3
            )
          )

          ###______________________________________________________________###
          ### 5) DBH estimation is done using `TreeLS::tlsInventory`
          ###______________________________________________________________###
          ### Search through tree points and estimate DBH to return a data frame of results
            tree_inv_df <- TreeLS::tlsInventory(
              las = las_norm_tile
              # height layer (above ground) to estimate stem diameters, in point cloud units
              , dh = 1.37
              # height layer width, in point cloud units
              , dw = 0.5
              # height percentile to extract per tree (0-1). Use 1 for top height, i.e. the highest point.
              , hp = 1
              # parameterized shapeFit function, i.e. method to use for diameter estimation.
              , d_method = TreeLS::shapeFit(
                # either "circle" or "cylinder".
                shape = "circle"
                # optimization method for estimating the shape's parameters
                , algorithm = "ransac"
                # number of points selected on every RANSAC iteration.
                , n = 20
              )
            )
            # class(tree_inv_df)
            # tree_inv_df %>% dplyr::glimpse()
        if(nrow(tree_inv_df)>0){
          ###_______________________________________________________###
          ### 93) clean up the DBH stem data frame ###
          ###_______________________________________________________###
            # add details to table and convert to sf data
            tree_inv_df <- tree_inv_df %>%
              dplyr::mutate(
                Radius = as.numeric(Radius)
                , dbh_m = Radius*2
                , dbh_cm = dbh_m*100
                , basal_area_m2 = pi * (Radius)^2
                , basal_area_ft2 = basal_area_m2 * 10.764
                , treeID = paste0(X, "_", Y)
                , stem_x = X
                , stem_y = Y
              ) %>%
              sf::st_as_sf(coords = c("X", "Y"), crs = sf::st_crs(las_norm_tile)) %>%
              dplyr::select(
                treeID, H, stem_x, stem_y, Radius, Error
                , dbh_m, dbh_cm, basal_area_m2, basal_area_ft2
              ) %>%
              dplyr::rename(
                tree_height_m = H
                , radius_m = Radius
                , radius_error_m = Error
              )
            # tree_inv_df %>% dplyr::glimpse()

            ### Remove points outside the bounding box of the laz tile + 1m buffer
            tree_inv_df <- tree_inv_df %>%
              sf::st_crop(
                sf::st_bbox(las_norm_tile) %>%
                  sf::st_as_sfc() %>%
                  sf::st_buffer(1)
              )
          ### !!!!!!! Don't need to write these las files because they aren't used for anything ###
          # ###_______________________________________________________###
          # ### Set the classification codes of different point types ###
          # ###_______________________________________________________###
          #
          # ### Pull out the stem files
          # stem_points <- lidR::filter_poi(las_norm_tile, Stem == TRUE)
          # stem_points@data$Classification <- 4
          #
          # ### Pull out the ground points
          # ground <- lidR::filter_poi(las_norm_tile, Classification %in% c(2,9))
          #
          # ### Pull out the remaining points that arent ground
          # remaining_points <- lidR::filter_poi(las_norm_tile, Stem == FALSE & !(Classification %in% c(2,9)))
          # remaining_points@data$Classification <- 5
          #
          # ### Combine the newly classified data
          # las_reclassified <- rbind(stem_points, ground, remaining_points)
          # # str(las_reclassified)
          # # class(las_reclassified)
          # # plot(las_reclassified, color = "Classification")
          #
          # ###_______________________________________________________###
          # ### Write output to disk ###
          # ###_______________________________________________________###
          # ### Write the stem points to the disk
          # if(class(las_reclassified)=="LAS"){
          #   lidR::writeLAS(las_reclassified, paste0(normalizePath(outfolder), "/", las_name))
          # }

          ### Write stem polygons to the disk
          out_name <- tools::file_path_sans_ext(las_name)
          out_name <- paste0(normalizePath(outfolder), "/", out_name, ".gpkg")
          if(max(class(tree_inv_df)=="sf")==1){
            sf::st_write(obj = tree_inv_df, dsn = out_name, quiet = T, append = F)
          }
          return(T)
        }else{return(F)} # nrow(tree_inv_df)>0
      }else{return(F)} # treels_treemap() return is LAS
    }else{return(F)} # max_point_height >= min_height
  } # treels_pipeline

###___________________________________________###
### 3/3 chunk the                             ###
### normalized las if not already             ###
###___________________________________________###
  chunk_norm_las <- function(folder, outfolder) {
    ###################################
    # tile the normalized files to process with TreeLS
    ###################################
    # check if folder contains las files directly
    chk <- folder %>%
      tolower() %>%
      stringr::str_detect(".*\\.(laz|las)$")
    # file list
    if(max(chk)==1){
      ff <- stringr::str_subset(folder, ".*\\.(laz|las)$")
    }else{ # otherwise just search the folder
      ff <- list.files(normalizePath(folder), pattern = ".*\\.(laz|las)$", full.names = T)
    }
    # create spatial index files (.lax)
    flist_temp <- create_lax_for_tiles(las_file_list = ff)
    # read
    norm_ctg_temp <- lidR::readLAScatalog(flist_temp)
    # retile catalog
    lidR::opt_progress(norm_ctg_temp) <- F
    lidR::opt_output_files(norm_ctg_temp) <- paste0(normalizePath(outfolder),"/", "_{XLEFT}_{YBOTTOM}") # label outputs based on coordinates
    lidR::opt_chunk_buffer(norm_ctg_temp) <- 10
    lidR::opt_chunk_size(norm_ctg_temp) <- 100
    lidR::catalog_retile(norm_ctg_temp) # apply retile
    # create spatial index
    normalize_flist <- create_lax_for_tiles(
      las_file_list = list.files(normalizePath(outfolder), pattern = ".*\\.(laz|las)$", full.names = T)
    )
    # return
    return(normalize_flist)
  }
