#' @title Download Forest Type Groups of the Continental United States data
#' @param force Whether to overwrite existing data
#' @param savedir Optional directory to save data in a new location. Defaults to package contents.
#' @references
#' * [LANDFIRE Forest Canopy Bulk Density (CBD)](https://landfire.gov/fuel/cbd)
#' U.S. Department of Agriculture and U.S. Department of the Interior.
#'
#' @description
#' The LANDFIRE Forest Canopy Bulk Density (CBD) data is used to estimate individual tree crown biomass in kilograms
#' See [trees_biomass_landfire()]
#'
#' @examples
#'  \dontrun{
#'  get_landfire()
#'  }
#' @export
#'
get_landfire <- function(
  savedir = NULL
  , force = F
){
  # filename !!!!!!!!!!!!! landfire only
  ff_name <- "lc23_cbd_240.tif"
  # set up parameters to pass to get_url_data()
  my_eval_url <- "https://landfire.gov/data-downloads/US_240/LF2023_CBD_240_CONUS.zip"
  my_my_name <- "landfire"
  my_req_file_list <- c(ff_name)
  my_cleanup_zip <- T
  # set up to save csv to package contents with location of data
  # the package directory from get_url_data()
  pkg_dir <- pkg_dir()
  if(!dir.exists(pkg_dir)){
    dir.create(pkg_dir, showWarnings = FALSE)
  }
  my_savedir <- ifelse(
      purrr::is_empty( normalizePath(file.path(savedir)) )
      , pkg_dir
      , normalizePath(file.path(savedir))
    )

  # call get_url_data()
  get_ans <- get_url_data(
    eval_url = my_eval_url
    , my_name = my_my_name
    , savedir = my_savedir
    , req_file_list = my_req_file_list
    , force = force
    , cleanup_zip = my_cleanup_zip
  )

  # save the location to a csv file if successful download
  if(get_ans==T){
    # where was this written?
    fff <- file.path(my_savedir, my_my_name)
    # write a csv to package directory with location of data
    dplyr::tibble(location = fff) %>%
      write.csv(
        file.path(pkg_dir, "location_landfire.csv")
        , row.names = F
        , append = F
      )

    ##### !!!!!!!!! for landfire only...
    ##### !!!!!!!!! the original data is factor type
    ##### !!!!!!!!! the first time we download this data we'll clean it
    ##### !!!!!!!!! by converting the data to numeric and transforming the values
    # reclass_landfire_rast() defined in utils_rast_points.R
    fff_name <- file.path(fff, ff_name)
    safe_reclass_landfire_rast <- purrr::safely(reclass_landfire_rast)
    rcl <- safe_reclass_landfire_rast(rast = terra::rast(fff_name))
    if(is.null(rcl$error)){
      message("unpacking LANDFIRE data...this might take a while (24-98 mins.)")
      terra::writeRaster(rcl$result, filename = fff_name, overwrite = T)
      message("LANDFIRE raster successfully downloaded and unpacked")
    }else{
      warning("LANDFIRE raster successfully downloaded and but not fully unpacked...proceed with caution")
    }
  }
}
