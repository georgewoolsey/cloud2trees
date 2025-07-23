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
