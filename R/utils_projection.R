#' @title tools for working with projections, crs, epsg, and whatnot
#'
#' @param x sf. or any object of a class that can be read by sf::st_crs such as "LAS" class
#'
#' @return A data.frame with the units and scale factor for the vertical projection
#'
#' @keywords internal
#'
get_vertical_crs <- function(x) {
  xcrs <- sf::st_crs(x)
  if (is.na(xcrs)){
    # stop("No CRS defined...try setting the parameter `new_crs` if known")
    return(NA)
  }

  wkt <- sf::st_as_text(xcrs)

  if (!grepl("COMPD_CS", wkt)) {
    # Should just be a horizontal CRS - nothing
    message("no vertical crs found")
    return(NA)
  }

  # extract the VERT_CS component
  if(!stringr::str_detect(wkt,"VERT_CS")){
    message("no vertical crs found")
    return(NA)
  }

  # get the vertical component
  # might also be VERT_DATUM ???
  i <- regexpr("VERT_CS\\[", wkt)
  wkt <- base::substring(wkt, i)

  # find the units
  if(!stringr::str_detect(wkt,"UNIT")){
    message("no vertical crs found")
    return(NA)
  }

  i <- regexpr("UNIT\\[", wkt)
  wkt <- base::substring(wkt, i) %>% stringr::str_remove("UNIT\\[")

  # work with the unit part
  unit_part <- stringr::str_extract(wkt, "[^\\]]+")
  unit_df <-
    read.csv(text = unit_part, header = F, colClasses = "character") %>%
    dplyr::mutate(
      units = stringr::str_extract(
        tolower(V1)
        , paste(c("foot","feet","meter","metere","metre"),collapse = "|")
      )
      , scale_factor = readr::parse_number(V2)
    ) %>%
    dplyr::select(units,scale_factor) %>%
    dplyr::slice(1)

  if(nrow(unit_df)==1){
    return(unit_df)
  }else{
    return(NA)
  }

}

###___________________________________________###
# Retrieve the horizontal component of a compound CRS.
# The object x can be an 'sf' package 'crs' object or any
# spatial object from which a CRS can be queried using the
# sf::st_crs function.
###___________________________________________###
get_horizontal_crs <- function(x) {
  xcrs <- sf::st_crs(x)
  if (is.na(xcrs)){
    # stop("No CRS defined...try setting the parameter `new_crs` if known")
    return(NA)
  }

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

###___________________________________________###
# transform raw las x,y points and return a data.frame
# of x,y in new projection which can be combined with
# a Z column for full definition of LAS
# https://gis.stackexchange.com/questions/375735/clean-way-to-convert-from-us-feet-to-meters-with-lidr-lidar-point-cloud
# https://github.com/r-lidar/lidR/issues/372
###___________________________________________###
st_transform_las_xy <- function(las, new_epsg_code) {
  stopifnot(inherits(las, "LAS"))
  # get epsg number
  if(inherits(new_epsg_code,"character")){
    e <- readr::parse_number(new_epsg_code)
  }else{
    e <- new_epsg_code
  }
  if(
    is.na(e) || is.null(e) ||
    identical(e,character(0)) ||
    !inherits(e,"numeric")
  ){
    stop("could not detect numeric epsg code")
  }
  # do the work and return data.frame with X,Y that can be converted to las after adding Z
  if(lidR::is.empty(las)){return(NULL)}
  return_df <- las@data %>%
    dplyr::select(X,Y) %>%
    sf::st_as_sf(
      coords = c("X", "Y")
      , crs = lidR::st_crs(las)
      , remove = F
    ) %>%
    dplyr::rename(orig_x = X, orig_y = Y) %>%
    sf::st_transform(crs = e) %>%
    dplyr::mutate(
      X = sf::st_coordinates(.)[,1]
      , Y = sf::st_coordinates(.)[,2]
    ) %>%
    dplyr::select(X,Y) %>%
    sf::st_drop_geometry()
  # return
  return(return_df)
}

###___________________________________________###
# transform raw las z points and return a data.frame
# of z scaled which can be combined with
# a x,y column for full definition of LAS
# https://gis.stackexchange.com/questions/375735/clean-way-to-convert-from-us-feet-to-meters-with-lidr-lidar-point-cloud
# https://github.com/r-lidar/lidR/issues/372
###___________________________________________###
st_transform_las_z <- function(
  las
  , scale_factor # this gets multiplied by current z
){
  stopifnot(inherits(las, "LAS"))
  # get epsg number
  if(inherits(scale_factor,"character")){
    f <- readr::parse_number(scale_factor)
  }else{
    f <- scale_factor
  }
  if(
    is.na(f) || is.null(f) ||
    identical(f,character(0)) ||
    !inherits(f,"numeric")
  ){
    stop("could not detect numeric scale factor")
  }

  # do the work and return data.frame with z that can be converted to las after adding x,y
  if(lidR::is.empty(las)){return(NULL)}
  return_df <- las@data %>%
    dplyr::select(Z) %>%
    dplyr::mutate(
      Z = Z*f
    )
  # return
  return(return_df)
}
