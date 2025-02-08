#' @title Tile raw `.las`|`.laz` files to work with smaller chunks
#'
#' @param folder string. The path of a folder containing a set of las/laz files. Can also be a vector of file paths.
#' @param outfolder string. The path of a folder to write the tiled las files to.
#' @param accuracy_level numeric. Choose processing accuracy.
#'      accuracy_level = 1 uses DTM to height normalize the points
#'      accuracy_level = 2 uses triangulation with high point density (20 pts/m2) to height normalize the points
#'      accuracy_level = 3 uses triangulation with very high point density (100 pts/m2) to height normalize the points
#' @param max_ctg_pts numeric. Max number of points to process at one time. Setting this number higher will possibly reduce run times but increase the chance of running out of memory and vice versa.
#' @param max_area_m2 numeric. Max area to process at one time. See `max_ctg_pts` parameter, this one is less important as never experienced memory issues with large areas (just lots of points)
#' @param transform logical. should the las/laz files be transformed? If set to `TRUE` the parameters `new_crs` must be defined.
#' @param new_crs string. crs to change to as an epsg numerical code
#' @param old_crs string. crs to change from as an epsg numerical code
#'
#' @references
#' [https://r-lidar.github.io/lidRbook/norm.html](https://r-lidar.github.io/lidRbook/norm.html)
#' [https://github.com/r-lidar/lasR/issues/18#issuecomment-2027818414](https://github.com/r-lidar/lasR/issues/18#issuecomment-2027818414)
#'
#' @description
#' Function to tile raw `.las`|`.laz` files to work with smaller chunks based on point density and coverage area
#'
#' @return A list of 1) `process_data` an sf object; 2) `is_chunked_grid` indicator if chunks were created; 3) `plt` a ggplot object
#'
#' @examples
#'  \dontrun{
#'  f <- "../lasdata"
#'  chunk_las_catalog(folder = f, outfolder = getwd())
#'  }
#' @export
#'
chunk_las_catalog <- function(
  folder
  , outfolder = getwd()
  , accuracy_level = 2
  , max_ctg_pts = 70e6
  , max_area_m2 = 90e6
  , transform = FALSE
  , new_crs = NA
  , old_crs = NA
){
  # accuracy level
  #### This maximum filters the ground points to perform Delaunay triangulation
  # see: https://github.com/r-lidar/lasR/issues/18#issuecomment-2027818414
  max_pts_m2 <- dplyr::case_when(
    as.numeric(accuracy_level) <= 2 ~ 20
    , as.numeric(accuracy_level) == 3 ~ 100
    , T ~ 20
  )
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
  ### point to input las files as a lidR LAScatalog (reads the header of all the LAS files of a given folder)
    las_ctg <- lidR::readLAScatalog(flist_temp)

  ###______________________________###
  # crs checks and transformations
  ###______________________________###
    # transform if user requested
    if(
      transform == TRUE &
      !is.na(new_crs) &
      !is.null(new_crs) &
      length(as.character(new_crs))>0
    ){
      message("reprojecting raw las data as requested. may lead to result inaccuracies!")
      # map over files
        flist_temp <- las_ctg@data$filename %>%
          purrr::map(
            reproject_las
            , new_crs = new_crs
            , old_crs = old_crs
            , outdir = normalizePath(outfolder)
          ) %>%
          c() %>%
          unlist()
      # create spatial index files (.lax)
        flist_temp <- create_lax_for_tiles(flist_temp)
      ### point to input las files as a lidR LAScatalog
        las_ctg <- lidR::readLAScatalog(flist_temp)
    }

    # pull crs for using in write operations
    # crs_list_temp = las_ctg@data$CRS
    crs_list_temp <- sf::st_crs(las_ctg)$epsg

    # handle missing epsg with user defined parameter
    if(is.na(crs_list_temp) & !is.na(new_crs) & !is.null(new_crs)){
      crs_list_temp <- new_crs
      sf::st_crs(las_ctg) <- paste0("EPSG:", new_crs)
    }else if(is.na(crs_list_temp)){
      # try to pull the epsg another way if still NA
      n_crs <- get_horizontal_crs(las_ctg)
      if(!is.na(n_crs)){
        las_ctg <- las_ctg %>% sf::st_set_crs(n_crs)
        crs_list_temp <- sf::st_crs(las_ctg)$epsg
      }else{
        stop("No CRS defined...try setting the parameter `new_crs` if known")
      }
    }
    # set crs for use in project
    if(length(unique(crs_list_temp))>1){
      stop("The raw las files have multiple CRS settings. Confine las files in `folder` to files with same CRS or re-generate las files with same projection.")
    }else{
      proj_crs <- paste0("EPSG:",unique(crs_list_temp))
    }

  ### is this ctg huge or what?
    ctg_pts_so_many <- sum(las_ctg@data$Number.of.point.records) > max_ctg_pts

  ########################
  # check for tile overlaps and so many points in whole catalog
  # ...determines tiling and processing of grid subsets
  ########################
    ctg_chunk_data <-
      las_ctg@data %>%
      dplyr::mutate(
        area_m2 = sf::st_area(geometry) %>% as.numeric()
        , pts = Number.of.point.records
      ) %>%
      sf::st_drop_geometry() %>%
      dplyr::select(filename, area_m2, pts) %>%
      dplyr::ungroup() %>%
      dplyr::summarise(
        dplyr::across(.cols = c(area_m2, pts), .fns = sum)
      ) %>%
      # attach overall area
      dplyr::bind_cols(
        las_ctg@data$geometry %>%
          sf::st_union() %>%
          sf::st_as_sf() %>%
          sf::st_set_geometry("geometry") %>%
          dplyr::mutate(
            total_area_m2 = sf::st_area(geometry) %>% as.numeric()
          ) %>%
          sf::st_drop_geometry()
      ) %>%
      dplyr::mutate(
        # so many points chunk size
        # uneven distribution of points in ctg => some chunks bigger than max_ctg_pts
        ctg_pt_factor = pts/(max_ctg_pts*0.2) #...so make much smaller
        , chunk_max_ctg_pts = dplyr::case_when(
            # no resize
            ctg_pt_factor <= 1 ~ 0
            # yes resize
            , ctg_pt_factor > 1 ~ round(sqrt(area_m2/ctg_pt_factor), digits = -1) # round to nearest 10
          )
        , buffer_chunk_max_ctg_pts = ifelse(round(chunk_max_ctg_pts*0.05,digits=-1)<10,10,round(chunk_max_ctg_pts*0.05,digits=-1))
        # overlap chunk size
        , pct_overlap = (area_m2-total_area_m2)/total_area_m2
        , area_factor = total_area_m2/max_area_m2 ## if want to split in to area-based chunks, define max_area_m2 parameter
        , pt_factor = pts/max_ctg_pts
        , chunk_overlap = dplyr::case_when(
            # no resize
            area_factor <= 1 & pt_factor <= 1 ~ 0
            # yes resize
            , area_factor > 1 & pt_factor > 1 ~ min(
              round(sqrt(area_m2/area_factor), digits = -1) # round to nearest 10
              , round(sqrt(area_m2/pt_factor), digits = -1) # round to nearest 10
            )
            , area_factor > 1 ~ round(sqrt(area_m2/area_factor), digits = -1) # round to nearest 10
            , pt_factor > 1 ~ round(sqrt(area_m2/pt_factor), digits = -1) # round to nearest 10
          )
        , buffer_overlap = ifelse(round(chunk_overlap*0.05,digits=-1)<10,10,round(chunk_overlap*0.05,digits=-1))
      )
    # ctg_chunk_data
  ############################################################
  ############################################################
  # use ctg_chunk_data to determine chunking
  ############################################################
  ############################################################
    # plot(las_ctg)
    ########################
    # split into grid subsets
    # this is done for ctg's with sooo many points
    # because lasR pipeline crashes memory and only way around it is
    # to break out las into separarte chunks (requires reading/writing much data ;()
    # ....
    # objective here is to break the ctg into lasR manageable subsets: each subset goes through lasR pipeline individually
    # then bring the DTM and CHM results together via terra::sprc and mosaic, e.g.:
      # # read
      # rast_list = list.files("../data/", pattern = ".*\\.(tif|tiff)$", full.names = T) %>% purrr::map(function(x){terra::rast(x)})
      # # mosaic
      # rast_mosaic = terra::sprc(rast_list) %>% terra::mosaic(fun = "max")
    ########################

    ########################
    # first, set the lidR select option
    ########################
    # # # from the lidR documentation
        # # # the select argument specifies the data that will actually be loaded.
        # #   For example, ’xyzia’ means that the x, y, and z coordinates, the intensity and the scan angle will be loaded.
        # #   The supported entries are t - gpstime, a - scan angle, i - intensity, n - number of returns, r - return number
        # #   , c - classification, s - synthetic flag, k - keypoint flag, w - withheld flag, o - overlap flag (format 6+)
        # #   , u - user data, p - point source ID, e - edge of flight line flag, d - direction of scan flag
        # #   , R - red channel of RGB color, G - green channel of RGB color, B - blue channel of RGB color
        # #   , N - near-infrared channel, C - scanner channel (format 6+), W - Full waveform.
        # #   Also numbers from 1 to 9 for the extra bytes data numbers 1 to 9. 0 enables all extra bytes to
        # #     be loaded and ’*’ is the wildcard that enables everything to be loaded from the LAS file.
        lidr_select <- "xyzainrcRGBNC"

    if(
      ctg_chunk_data$chunk_max_ctg_pts[1] > 0
      & ctg_pts_so_many == T
    ){
      # if need to retile
      # retile catalog
      lidR::opt_progress(las_ctg) <- F
      lidR::opt_output_files(las_ctg) <- paste0(normalizePath(outfolder),"/", "_{XLEFT}_{YBOTTOM}") # label outputs based on coordinates
      lidR::opt_filter(las_ctg) <- "-drop_duplicates"

      # # https://gis.stackexchange.com/questions/378882/how-can-i-make-las-data-output-from-an-rpas-survey-processed-in-agisoft-metashap
      lidR::opt_select(las_ctg) <- lidr_select
      # # lidR::opt_select(las_ctg) <- "-x -y -z -i -n -r -c"
      # # lidR::opt_select() <- "-u -i -w"

      # buffering here because these grid subsets will be processed independently
      lidR::opt_chunk_buffer(las_ctg) <- 10
      lidR::opt_chunk_size(las_ctg) <- ctg_chunk_data$chunk_max_ctg_pts[1]
      # reset las_ctg
      lidR::catalog_retile(las_ctg) # apply retile
      lidR::opt_progress(las_ctg) <- T
      # create spatial index
      flist_temp <- create_lax_for_tiles(
        las_file_list = list.files(normalizePath(outfolder), pattern = ".*\\.(laz|las)$", full.names = T)
      )
      # switch for processing grid subsets
      is_chunked_grid <- T

    }else if( # retile whole catalog if high overlap
      ctg_chunk_data$chunk_overlap[1] > 0
      & ctg_chunk_data$pct_overlap[1] > 0.1
    ){
      # if need to retile
      # retile catalog
      lidR::opt_progress(las_ctg) <- F
      lidR::opt_output_files(las_ctg) <- paste0(normalizePath(outfolder),"/", "_{XLEFT}_{YBOTTOM}") # label outputs based on coordinates
      lidR::opt_filter(las_ctg) <- "-drop_duplicates"
      # # # from the lidR documentation
        # # # the select argument specifies the data that will actually be loaded.
        # #   For example, ’xyzia’ means that the x, y, and z coordinates, the intensity and the scan angle will be loaded.
        # #   The supported entries are t - gpstime, a - scan angle, i - intensity, n - number of returns, r - return number
        # #   , c - classification, s - synthetic flag, k - keypoint flag, w - withheld flag, o - overlap flag (format 6+)
        # #   , u - user data, p - point source ID, e - edge of flight line flag, d - direction of scan flag
        # #   , R - red channel of RGB color, G - green channel of RGB color, B - blue channel of RGB color
        # #   , N - near-infrared channel, C - scanner channel (format 6+), W - Full waveform.
        # #   Also numbers from 1 to 9 for the extra bytes data numbers 1 to 9. 0 enables all extra bytes to
        # #     be loaded and ’*’ is the wildcard that enables everything to be loaded from the LAS file.

      # # https://gis.stackexchange.com/questions/378882/how-can-i-make-las-data-output-from-an-rpas-survey-processed-in-agisoft-metashap      # https://gis.stackexchange.com/questions/378882/how-can-i-make-las-data-output-from-an-rpas-survey-processed-in-agisoft-metashap
      lidR::opt_select(las_ctg) <- lidr_select
      # # lidR::opt_select(las_ctg) <- "-x -y -z -i -n -r -c"
      # # lidR::opt_select() <- "-u -i -w"

      # buffering is handled by lasR::exec so no need to buffer here
      # not buffering here because these grid subsets will be processed altogether with buffer set for lasR::exec
      lidR::opt_chunk_buffer(las_ctg) <- 0
      lidR::opt_chunk_size(las_ctg) <- ctg_chunk_data$chunk[1]
      # reset las_ctg
      lidR::catalog_retile(las_ctg) # apply retile
      lidR::opt_progress(las_ctg) <- T
      # create spatial index
      flist_temp <- create_lax_for_tiles(
        las_file_list = list.files(normalizePath(outfolder), pattern = ".*\\.(laz|las)$", full.names = T)
      )
      # switch for processing grid subsets
      is_chunked_grid <- F

    }else{
      is_chunked_grid <- F

    }
      # flist_temp
      # is_chunked_grid
      # plot(lidR::readLAScatalog(flist_temp))

      # lidR::readLAScatalog(flist_temp)@data %>%
      #   dplyr::select(filename, Number.of.point.records) %>%
      #   sf::st_drop_geometry()

  ########################
  # data on how the chunks are processed for writing
  # get pct filter for creating dtm and normalized point clouds using triangulation
  ########################
    process_data <-
      lidR::readLAScatalog(flist_temp)@data %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        processing_grid = dplyr::case_when(
          is_chunked_grid == T ~ dplyr::row_number()
          , T ~ 1
        )
        , area_m2 = sf::st_area(geometry) %>% as.numeric()
        , pts = Number.of.point.records
        , pts_m2 = pts/area_m2
      ) %>%
      dplyr::select(processing_grid, filename, area_m2, pts, pts_m2) %>%
      # get rid of tiny edge files with very few points which cannot be utilized for delauny triangulation
        # !!!!!!!!!!!!!!! upgrade to remove tiles completely spatially covered by other tiles
      dplyr::filter(
        ( pts>max_ctg_pts*0.0001 & dplyr::n() > 1 ) ## only remove small chunks if more than one
        | (dplyr::n()==1)
      ) %>%
      # sometimes there are super small chunks created which can lead to 0 points after filtering for performing calcs
      dplyr::mutate(
        processing_grid = dplyr::case_when(
          pts < max_ctg_pts*0.0002 & !is.na(dplyr::lag(pts)) & !is.na(dplyr::lead(pts)) &
            dplyr::lag(pts) < dplyr::lead(pts) ~ dplyr::lag(processing_grid)
          , pts < max_ctg_pts*0.0002 & !is.na(dplyr::lag(pts)) & !is.na(dplyr::lead(pts)) &
            dplyr::lag(pts) > dplyr::lead(pts) ~ dplyr::lead(processing_grid)
          , pts < max_ctg_pts*0.0002 & !is.na(dplyr::lag(pts)) ~ dplyr::lag(processing_grid)
          , pts < max_ctg_pts*0.0002 & !is.na(dplyr::lead(pts)) ~ dplyr::lead(processing_grid)
          , T ~ processing_grid
        )
      ) %>%
      dplyr::group_by(processing_grid) %>%
      dplyr::mutate(
        pts_m2 = pts/area_m2
      # calculate factor to reduce by for triangulation
        , pts_m2_factor_optimal = dplyr::case_when(
            max_pts_m2/pts_m2 >= 1 ~ 1
            , TRUE ~ max_pts_m2/pts_m2
        )
        # area based weighted mean
          # pass this to filter_for_dtm
        , pts_m2_factor_ctg = max(stats::weighted.mean(
            x = pts_m2_factor_optimal
            , w = area_m2/sum(area_m2)
          ) %>%
          round(digits = 2),0.01)
        , filtered_pts_m2 = pts_m2*pts_m2_factor_ctg
        , processing_grid_tot_pts = sum(pts)
        , normalization_accuracy = accuracy_level
        , proj_crs = proj_crs
      ) %>%
      dplyr::ungroup()

    # plot
    plt <- process_data %>%
      dplyr::group_by(processing_grid) %>%
      dplyr::mutate(lab = pts == max(pts)) %>%
      ggplot2::ggplot() +
        ggplot2::geom_sf(ggplot2::aes(fill = as.factor(processing_grid)), alpha = 0.6) +
        ggplot2::geom_sf_text(ggplot2::aes(label = ifelse(lab,as.factor(processing_grid), ""))) +
        ggplot2::scale_fill_manual(
          values =
            viridis::turbo(n=length(process_data$processing_grid %>% unique())) %>%
            sample()
        ) +
        ggplot2::labs(x="",y="",title="processing_grid") +
        ggplot2::theme_light() +
        ggplot2::theme(legend.position = "none")
    # # count
    # process_data %>%
    #   dplyr::distinct(processing_grid, processing_grid_tot_pts)

  # # save processing attributes
  #   sf::st_write(
  #     process_data
  #     , paste0(config$delivery_dir, "/process_data_las_tiling.gpkg")
  #     , quiet = TRUE, append = FALSE
  #   )

  # return
    return(list(
      process_data = process_data
      , is_chunked_grid = is_chunked_grid
      , plt = plt
      , las_ctg = las_ctg
    ))
}
###_____________________________________________________###
### Intermediate functions ###
###_____________________________________________________###

###___________________________________________###
# Retrieve the horizontal component of a compound CRS.
# The object x can be an 'sf' package 'crs' object or any
# spatial object from which a CRS can be queried using the
# sf::st_crs function.
###___________________________________________###
get_horizontal_crs <- function(x) {
  xcrs <- sf::st_crs(x)
  if (is.na(xcrs)) stop("No CRS defined...try setting the parameter `new_crs` if known")

  wkt <- sf::st_as_text(xcrs)

  if (!grepl("COMPD_CS", wkt)) {
    # Should just be a horizontal CRS - simply return it
    xcrs
  } else {
    # Extract the horizontal component
    i <- regexpr("PROJCS\\[", wkt)
    wkt <- base::substring(wkt, i)

    # Match square brackets to discard any trailing
    # component (e.g. the vertical CRS)
    wkt_chars <- base::strsplit(wkt, "")[[1]]
    level <- 1
    k <- base::match("[", wkt_chars)
    while (level > 0) {
      k <- k + 1
      if (wkt_chars[k] == '[') {
        level <- level + 1
      } else if (wkt_chars[k] == ']') {
        level <- level - 1
      }
    }

    wkt <- base::substring(wkt, 1, k)
    # return
    return(sf::st_crs(wkt))
  }
}
