#' @title Individual Tree Detection (ITD) tuning
#'
#' @description
#' `itd_tuning()` is used to visually assess tree crown delineation results
#' from different window size functions used for the detection of individual trees.
#' The `cloud2trees` package performs individual tree detection using [lidR::locate_trees()] with the [lidR::lmf()] algorithm.
#' The local maximum filter algorithm allows for a constant window size or a variable window size defined by a function.
#' See the `lidR` [package book](https://r-lidar.github.io/lidRbook/itd.html) for excellent detail on ITD and defining window size.
#'
#' `itd_tuning()` allows users to test different window size functions on a sample of data to determine which
#' function is most suitable for the area being analyzed. The preferred function can then be used in the `ws`
#' parameter in [raster2trees()] and/or [cloud2trees()].
#'
#'
#' @param input_las_dir character. directory where .las|.laz point cloud data exists...program will search all sub-directories for all .las|.laz files and process them as one
#' @param n_samples numeric. The number of sample plots of 0.1 ha on which to test the window functions. The maximum is 5.
#' The center of the point cloud data coverage will always be the first plot sampled so long as points exist in the central 0.1 ha.
#' @param ws_fn_list list. A function or a named list of functions. Leave as NULL to test default exponential (concave up), linear, and logarithmic (concave down) functions.
#' If providing a custom function, it must always return a numeric value >0 (see examples).
#' @param min_height numeric. Set the minimum height (m) for individual tree detection
#' @param chm_res_m numeric. The desired resolution of the CHM produced in meters.
#'
#' @references
#' [https://r-lidar.github.io/lidRbook/itd.html](https://r-lidar.github.io/lidRbook/itd.html)
#'
#' @return Returns a list with: 1) "plot_samples" is a plot of the sample canopy height model (CHM) and extracted tree crowns for each window size tested;
#' and 2) "ws_fn_list" is a list of the window size functions tested which can be used to pass the desired function to the `ws`
#' parameter in [raster2trees()] and/or [cloud2trees()].
#'
#' @examples
#'  \dontrun{
#'   # do it
#'   library(tidyverse)
#'   # test las file but this could also be a directory path with >1 .las|.laz files
#'   i <- system.file(package = "lidR", "extdata", "MixedConifer.laz")
#'   ####################################################
#'   # check the default itd_tuning() window functions
#'   ####################################################
#'    # run it with defaults
#'    itd_tuning_ans <- itd_tuning(input_las_dir = i)
#'    # what's in it?
#'    names(itd_tuning_ans)
#'    # look at the tuning plot
#'    itd_tuning_ans$plot_samples
#'    # the "exp_fn" looks pretty good, let's store it
#'    best_default <- itd_tuning_ans$ws_fn_list$exp_fn
#'    # we can see what this function looks like for window size
#'    ggplot2::ggplot() +
#'      ggplot2::geom_function(fun = best_default) +
#'      ggplot2::xlim(-5,60) +
#'      ggplot2::labs(x = "heights", y = "ws", color = "")
#'    # pass our best function to the cloud2trees() to process the full point cloud coverage
#'    cloud2trees_ans <- cloud2trees(output_dir = tempdir(), input_las_dir = i, ws = best_default)
#'    # the same plot as the the tuning plot with tree crowns overlaid on CHM
#'    ggplot2::ggplot() +
#'      ggplot2::geom_tile(
#'        data = cloud2trees_ans$chm_rast %>%
#'          terra::as.data.frame(xy=T) %>%
#'          dplyr::rename(f=3)
#'        , mapping = ggplot2::aes(x = x, y = y, fill = f)
#'        , na.rm = T
#'      ) +
#'      ggplot2::scale_fill_viridis_c(
#'        option = "plasma"
#'        , breaks = scales::breaks_extended(n=10)
#'      ) +
#'      ggplot2::geom_sf(
#'        data = cloud2trees_ans$crowns_sf
#'        , fill = NA, color = "gray33", lwd = 1
#'      ) +
#'      ggplot2::scale_x_continuous(expand = c(0, 0)) +
#'      ggplot2::scale_y_continuous(expand = c(0, 0)) +
#'      ggplot2::labs(x = "", y = "", fill = "CHM (m)") +
#'      ggplot2::theme_light() +
#'      ggplot2::theme(axis.text = ggplot2::element_blank())
#'   ####################################################
#'   # let's test some custom window functions
#'   ####################################################
#'     # a constant window size has to be defined as:
#'      ## x*0 + constant
#'      my_constant <- function(x){(x * 0) + 3} ## will always return 3
#'     # a custom linear function
#'      my_linear <- function(x) {(x * 0.1) + 3}
#'     # run it with custom functions
#'      itd_tuning_ans2 <- itd_tuning(
#'        input_las_dir = i
#'        , ws_fn_list = list(
#'           my_constant=my_constant
#'           , my_linear=my_linear
#'           , best_default=best_default # the best from our first test
#'         )
#'        , n_samples = 2
#'       )
#'     # look at the tuning plot
#'      itd_tuning_ans2$plot_samples
#'     # we can see what our custom "my_linear" function looks like
#'      ggplot2::ggplot() +
#'        ggplot2::geom_function(fun = itd_tuning_ans2$ws_fn_list$my_linear) +
#'        ggplot2::xlim(-5,60) +
#'        ggplot2::labs(x = "heights", y = "ws", color = "")
#'  }
#' @export
#'
itd_tuning <- function(
  input_las_dir
  , n_samples = 3
  , ws_fn_list = NULL
  , min_height = 2
  , chm_res_m = 0.25
){
#####################################################
# check n_samples
#####################################################
  n_samples <- n_samples %>%
    as.numeric() %>%
    dplyr::coalesce(3) %>%
    min(5)
#####################################################
# check the las and read as a LASCatalog if exists
#####################################################
  las <- check_las_data(input_las_dir)
  # set the lascatalog options
  if(inherits(las, "LAScatalog")){
    lidR::opt_progress(las) <- F
    lidR::opt_filter(las) <- "-drop_duplicates" ## class 2 = ground; 9 = water; 18 = noise
    lidR::opt_select(las) <- "xyz" # 0 enables all extra bytes to be loaded...possibly treeID
    lidR::opt_output_files(las) <- paste0(tempdir(), "/{*}_treed")
  }
#####################################################
# bbox based on the points
#####################################################
  roi_poly <- las %>% lidR::st_bbox() %>% sf::st_as_sfc()
  # roi_poly %>% sf::st_area()
  # make the roi 5% smaller to sample from so we don't get sample points on the edge
  roi_poly <- roi_poly %>%
    sf::st_buffer(
      roi_poly %>%
      sf::st_area() %>%
      as.numeric() %>%
      sqrt() %>%
      `*`(-0.0125) %>%
      round()
    )
  # roi_poly %>% sf::st_area()
#####################################################
# data frame of sample points that fall within roi
#####################################################
  sample_pts <-
    # the centroid will always be the first sample point
    roi_poly %>%
      sf::st_centroid() %>%
      sf::st_as_sf() %>%
      dplyr::bind_rows(
        sf::st_sample(
          roi_poly
          , size = max(5,n_samples*2)
          , type = "random"
        ) %>%
        sf::st_as_sf()
      ) %>%
      dplyr::mutate(rnk = dplyr::row_number())

  # ggplot() +
  #   geom_sf(data = roi_poly, fill = NA) +
  #   geom_sf(data = sample_pts, aes(color = rnk)) +
  #   theme_void()

#####################################################
# map sample points over sample_las_point_extract_trees()
#####################################################
quiet_sample_las_point_extract_trees <- purrr::quietly(sample_las_point_extract_trees)
big_ans_list <- 1:nrow(sample_pts) %>%
  purrr::map(\(x) quiet_sample_las_point_extract_trees(
      las = las
      , sample_point = sample_pts[x,]
      , ws_fn_list = ws_fn_list
      , min_height = min_height
      , chm_res_m = chm_res_m
    )
  )
# get the result
big_ans_list <- big_ans_list %>% purrr::map(purrr::pluck("result"))
if(is.null(big_ans_list) || dplyr::coalesce(length(big_ans_list), 0)<1){
  stop("failed to extract trees in any sample. try increasing `n_samples` or using different window functions.")
}
# str(big_ans_list)
#####################################################
# get the first n answers that are complete
#####################################################
  has_ans <- big_ans_list %>%
    purrr::map(purrr::pluck("crowns")) %>%
    purrr::map(function(x) ifelse(is.null(x),F,nrow(x)>0)) %>%
    unlist()
  # big_ans_list[] %>% length()
  # big_ans_list[has_ans] %>% length()

  # fitler or filter big list
  small_ans_list <- big_ans_list[has_ans]
  if(is.null(small_ans_list) || dplyr::coalesce(length(small_ans_list), 0)<1){
    stop("failed to extract trees in any sample. try increasing `n_samples` or using different window functions.")
  }
  # get the desired number
  the_ans_list <- small_ans_list[1:min(n_samples,length(small_ans_list))]
  if(is.null(the_ans_list) || dplyr::coalesce(length(the_ans_list), 0)<1){
    stop("failed to extract trees in any sample. try increasing `n_samples` or using different window functions.")
  }
# pull out the plots and combine them
  plt_list <- the_ans_list %>% purrr::map(purrr::pluck("plt"))
  if(any(is.null(plt_list))){
    plt <- NULL
  }else{
    # the_ans_list %>% purrr::map(names)
    plt <- patchwork::wrap_plots(
        plt_list
        , nrow = length(the_ans_list)
        , ncol = 1
      ) +
      patchwork::plot_annotation(tag_levels = "1", tag_prefix = "sample #")
  }

# pull out the ws fn...will be the same for each list so just need to get one
  if(any(is.null(plt_list))){
    ws_fn_list <- NULL
  }else{
    ws_fn_list <- the_ans_list[1] %>%
      purrr::map(purrr::pluck("ws_fn_list")) %>%
      purrr::flatten()
  }

  return(list(
    plot_samples = plt
    , ws_fn_list = ws_fn_list
  ))
}

#################################################################
#################################################################
## intermediate function (actually the workhorse) to:
  ## clip the las based on a buffered point
  ## run cloud2trees::cloud2raster() to get chm
  ## run cloud2trees::raster2trees() over a list of ws functions
  ## return chm, crown polys unique by ws fn, and ws functions used
#################################################################
#################################################################
sample_las_point_extract_trees <- function(
  las # a LASCatalog, maybe a LAS
  , sample_point = NULL
  , plot_area_m2 = 1000
  , ws_fn_list = NULL
  ### ws_fn_list = a function or a named list of functions
  ### leave as NULL to test default exponential and linear functions
  ### if providing a custom function, it must always return a numeric value >0.
  , min_height = 2
  , chm_res_m = 0.25
) {
  #check las
  if(
    !inherits(las, "LAScatalog") && !inherits(las,"LAS")
  ){stop("las must be LAS or LAScatalog")}
  #####################################################
  # make an empty return object if steps don't work
  #####################################################
  empty_return <- list(
      crowns = NULL
      , chm_rast = NULL
      , plt = NULL
      , ws_fn_list = NULL
    )
  #####################################################
  # check ws functions
  #####################################################
  # default function list
    # ws_fn_list = NULL ### !!! testing only !!! better take it out and don't forget
    def_fn_list <- list(
      lin_fn = function(x){
        y <- dplyr::case_when(
         is.na(x) ~ 0.001
         , x < 0 ~ 0.001
         , x < 2 ~ 1
         , x > 30 ~ 5
         , TRUE ~ 0.75 + (x * 0.14)
        )
        return(y)
      }
      , exp_fn = function(x) {
        y <- dplyr::case_when(
          is.na(x) ~ 1e-3 # requires non-null
          , x < 0 ~ 1e-3 # requires positive
          # , x < 2.5 ~ 1 # set lower bound
          # , x > 40 ~ 6.7  # set upper bound
          # , TRUE ~ exp( (0.0446*x) + (x^-0.555) ) # used gamma regression so exp the result
          #### NEW
          , x < 3.6 ~ 1.25 + x*0.15 # set lower bound
          , x > 32.5 ~ 5  # set upper bound
          , TRUE ~ exp( (0.0446*x) + (x^-0.555) ) # used gamma regression so exp the result
        )
        return(y)
      }
      , log_fn = function(x) {
        y <- dplyr::case_when(
          is.na(x) ~ 0.001
          , x < 0 ~ 0.001
          , x < 2 ~ 0.6
          , x > 26.5 ~ 5
          , TRUE ~ exp(-(3.5 * (1/x)) + (x^0.17)))
          # , TRUE ~ 2.6 * (-(exp(-0.08*(x-2)) - 1)) + 3 ## fn from lidR book
        return(y)
      }
    )
  # check ws_fn_list
    if(is.function(ws_fn_list)){
      ws_fn_list <- list(custom_fn = ws_fn_list)
    }else if(inherits(ws_fn_list,"list")){
      # filter list for functions
      ws_fn_list <- ws_fn_list %>%
        purrr::keep(is.function)
      # check whats left in the list
      if(
        any(is.na(ws_fn_list))
        || any(is.null(ws_fn_list))
        || length(ws_fn_list)<0
      ){
        ws_fn_list <- def_fn_list
      }
    }else{ # default
      ws_fn_list <- def_fn_list
    }

    # ggplot() +
    #   geom_function(fun = ws_fn_list$lin_fn, aes(color="lin_fn")) +
    #   geom_function(fun = ws_fn_list$exp_fn, aes(color="exp_fn")) +
    #   geom_function(fun = ws_fn_list$log_fn, aes(color="log_fn")) +
    #   xlim(-5,60) +
    #   labs(x = "Z", y = "ws", color = "")
  #####################################################
  # check sample point
  #####################################################
  # default sample point to centroid
  ## could also check that the sf objects are POINT...maybe another day
    if(!inherits(sample_point, "sf")){
      sample_point <- las %>%
        lidR::st_bbox() %>%
        sf::st_as_sfc() %>%
        sf::st_centroid() %>%
        sf::st_as_sf()
    }else{
      # just get the first element if someone passes >1
      sample_point <- sample_point[1,]
    }
  # plot_area_m2
    buffer <- plot_area_m2 %>%
      as.numeric() %>%
      dplyr::coalesce(0) %>%
      `/`(4) %>%
      sqrt() %>%
      max(
        sqrt(1000/4) ## numerator = desired plot size in m2
      )
  # plot
    plot_poly <- sample_point %>%
      sf::st_buffer(buffer, endCapStyle = "SQUARE") %>%
      dplyr::mutate(dummy=1)

  # ggplot() +
  #   geom_sf(
  #     data = las %>% lidR::st_bbox() %>% sf::st_as_sfc()
  #     , aes(color = "las bbox")
  #     , fill = NA
  #   ) +
  #   geom_sf(data = plot_poly, aes(color = "plot"), fill=NA) +
  #   geom_sf(data = sample_point, aes(color = "point"))

  #####################################################
  # clip and save file to pass to cloud2trees
  #####################################################

  # clip las
    # at a minimum we need to set the opt_output_files if ctg
    if(inherits(las, "LAScatalog")){
      lidR::opt_progress(las) <- F
      lidR::opt_filter(las) <- "-drop_duplicates"
      lidR::opt_select(las) <- "xyz"
      # "dummy" has to exist as a column name in the roi polygon data
      lidR::opt_output_files(las) <- paste0(tempfile(), "_{dummy}")
    }
    # clip it
    clip_las <- lidR::clip_roi(las, plot_poly)
    # if ctg, the clipped las is already written to disk and saved
    # if las was a LAS originally, clip_las will be a LAS and we need to writeLAS
    # of the las so that Las can retrieve the LAS
    if(inherits(clip_las, "LAScatalog")){
      # clip_las <- lidR::readLAS(files = clip_las@data$filename[1])
      las_fp <- clip_las@data$filename[1] ## should never have more than one element b/c only one polygon
    }else if(inherits(clip_las, "LAS")){
      # we need to write the las
      las_fp <- lidR::writeLAS(clip_las, file = file.path(tempfile(), "clip.las"))
    }else{return(empty_return)}

  #####################################################
  # cloud2trees that stuff
  #####################################################
    # get the raster
    safe_cloud2raster <- purrr::safely(cloud2raster)
    cloud2raster_ans <- safe_cloud2raster(
      input_las_dir = las_fp
      , output_dir = tempdir()
      , min_height = min_height # 1.37 ## dbh
      , chm_res_m = chm_res_m
    )
    # get the result
    if(is.null(cloud2raster_ans$error)){
      cloud2raster_ans <- cloud2raster_ans$result
    }else{
      return(empty_return)
    }
    # cloud2raster_ans$chm_rast %>% terra::plot()


    # map cloud2trees over each ws function
    safe_raster2trees <- purrr::safely(raster2trees)
    raster2trees_ans_list <-
      ws_fn_list %>%
      purrr::map(\(x)
        safe_raster2trees(
          chm_rast = cloud2raster_ans$chm_rast
          , outfolder = tempdir()
          , min_height = min_height # 1.37 ## dbh
          , ws = x
        )
      )
  # raster2trees_ans_list %>% str()

    # parse it to get it
    crowns <-
      raster2trees_ans_list %>%
      purrr::imap(parse_raster2trees_ans_list) %>%
      dplyr::bind_rows()

    if(is.null(crowns) || !inherits(crowns,"data.frame")){return(empty_return)}
  #####################################################
  # plot it
  #####################################################
    # crowns %>% dplyr::glimpse()
    # crowns %>% sf::st_drop_geometry() %>% dplyr::count(ws_fn)
    # plot the chm
    p_chm <- ggplot2::ggplot() +
      # adding the plot boundary should keep the plots at the same scale
      # for when we go to patchwork separate samples together
      ggplot2::geom_sf(data = plot_poly, fill = NA, lwd = 0) +
      ggplot2::geom_tile(
        data = cloud2raster_ans$chm_rast %>%
          terra::as.data.frame(xy=T) %>%
          dplyr::rename(f=3)
        , mapping = ggplot2::aes(x = x, y = y, fill = f)
        , na.rm = T
      ) +
      ggplot2::scale_fill_viridis_c(
        option = "plasma"
        , breaks = scales::breaks_extended(n=10)
      ) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(x = "", y = "", fill = "CHM (m)") +
      ggplot2::theme_light() +
      ggplot2::theme(axis.text = ggplot2::element_blank())
    # p_chm
    # add crowns by ws_fn to chm plot
    plt <- p_chm +
      ggplot2::geom_sf(
        data = crowns %>%
          dplyr::group_by(ws_fn) %>%
          dplyr::mutate(
            lab = paste0(
              ws_fn, "\n# trees = "
              , scales::comma(dplyr::n(), accuracy = 1)
            )
          )
        , fill = NA, color = "gray33", lwd = 0.8
      ) +
      ggplot2::facet_grid(cols = dplyr::vars(lab)) +
      ggplot2::theme(strip.text = ggplot2::element_text(color = "black"))

    # return it
    return(list(
      crowns = crowns
      , chm_rast = cloud2raster_ans$chm_rast
      , plt = plt
      , ws_fn_list = ws_fn_list
    ))

}
# # test tester
# xxx <- sample_las_point_extract_trees(
#   las = las
#   , sample_point = sample_pts[11,]
# )
# str(xxx)
################################################################
# this is a dumb one but...
# pull out the crown polygons from raster2trees_ans_list
# which is the answer from mapping a function list over the raster2trees()
################################################################
parse_raster2trees_ans_list <- function(x,y){
  nm <- y
  ans <- x %>% purrr::pluck("result")
  if(is.null(ans)){return(NULL)}
  # name the polygon file
  ans <- ans %>%
    dplyr::mutate(ws_fn = nm) %>%
    dplyr::relocate(ws_fn)
  return(ans)
}
