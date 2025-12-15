#' @title Check a for a directory with .laz|.las files,
#' the path of a file, or an object of class `LAS`|`LASCatalog`
#'
#' @description
#' Check a character for a directory with .laz|.las files,
#' the path of a file, or an object of class `LAS`|`LASCatalog`
#'
#' @param las character. a directory with .laz|.las files files,
#' the path of a single .laz|.las file, -or- an object of class `LAS`|`LASCatalog`
#'
#' @keywords internal
#'
check_las_data <- function(
  las
){
  las_msg <- paste0(
    "`las` must contain a directory with las files, the path of a .laz|.las file"
    , "\n, -or- an object of class `LAS`. Please update the `las` parameter."
  )
  # check the class
  if(inherits(las, "character")){
    if(
      !any(stringr::str_ends(las, ".*\\.(laz|las)$"))
      && length(las)==1
      && dir.exists(normalizePath(las))
    ){
      # try to read directory for las files
      fls <- list.files(normalizePath(las), pattern = ".*\\.(laz|las)$", full.names = TRUE)
      # stop it if no files
      if(length(fls)<1){stop(las_msg)}
      # read it
      nlas_ctg <- lidR::readLAScatalog(fls)
    }else if(any(stringr::str_ends(las, ".*\\.(laz|las)$"))){
      # read it
      nlas_ctg <- stringr::str_subset(las, pattern = ".*\\.(laz|las)$") %>%
        lidR::readLAScatalog()
    }else{
      stop(las_msg)
    }
  }else if(inherits(las, "LAS")){
    nlas_ctg <- las
  }else if(inherits(las, "LAScatalog")){
    nlas_ctg <- las
  }else{
    stop(las_msg)
  }

  return(nlas_ctg)
}

# make a function to check las ctg and get rid of empty files or non-valid geoms
check_las_ctg_empty <- function(las_ctg, pts = T, geoms = T) {
  if(!inherits(las_ctg, "LAScatalog")){stop("`las_ctg` must be of class LAScatalog")}
  # final flist
  if(pts && geoms){
    good_fls <- las_ctg@data %>%
      dplyr::ungroup() %>%
      sf::st_transform(5070) %>%
      sf::st_make_valid() %>%
      dplyr::filter(
        !sf::st_is_empty(.)
        & Number.of.point.records>0
      ) %>%
      dplyr::pull(filename) %>%
      unique()
  }else if(pts){
    good_fls <- las_ctg@data %>%
      dplyr::filter(
        Number.of.point.records>0
      ) %>%
      dplyr::pull(filename) %>%
      unique()
  }else if(geoms){
    good_fls <- las_ctg@data %>%
      dplyr::ungroup() %>%
      sf::st_transform(5070) %>%
      sf::st_make_valid() %>%
      dplyr::filter(
        !sf::st_is_empty(.)
      ) %>%
      dplyr::pull(filename) %>%
      unique()
  }else{
    good_fls <- unique(las_ctg@data$filename)
  }
  # check
  if(dplyr::coalesce(length(good_fls),0)==0){
    stop("no valid point cloud files found...all missing geometry or empty points")
  }

  # read las ctg
  ret_ctg <- check_las_data(good_fls)
  # huh
  if(!inherits(ret_ctg, "LAScatalog")){stop("could not detect .las|.laz files...all missing geometry or empty points")}
  return(ret_ctg)
}
