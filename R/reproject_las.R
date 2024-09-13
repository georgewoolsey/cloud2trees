#' function to reproject las
#' @description
#' Function to reproject las data using `lidR` but be careful! This is inefficient and potentially causes inaccuracies due to transformations (see reference).
#' @references
#' [https://gis.stackexchange.com/questions/371566/can-i-re-project-an-las-file-in-lidr](https://gis.stackexchange.com/questions/371566/can-i-re-project-an-las-file-in-lidr)
#'
#' @param filepath the full file path of the .las|.laz file
#' @param new_crs crs to change to as an epsg numerical code
#' @param old_crs crs to change from as an epsg numerical code
#' @param outdir where to save this re-projected .las file
#'
#' @return A character
#'
#' @keywords internal
#'
  reproject_las <- function(filepath, new_crs = NA, old_crs = NA, outdir = getwd()) {
    if(is.null(new_crs) | is.na(new_crs)){stop("the new_crs must be provided")}
    # read individual file
    las <- lidR::readLAS(filepath)
    if(
      (is.null(sf::st_crs(las)$epsg) | is.na(sf::st_crs(las)$epsg))
      & ( is.null(old_crs) | is.na(old_crs) | old_crs=="" )
    ){
      stop("the raw las file has missing CRS and cannot be transformed. try setting old_crs if known")
    }else{
      # transform if know old crs
      if(
        (is.null(sf::st_crs(las)$epsg) | is.na(sf::st_crs(las)$epsg))
        & !is.null(old_crs) & !is.na(old_crs) & length(as.character(old_crs))>0
      ){
        sf::st_crs(las) <- paste0("EPSG:", old_crs)
      }
      # get filename
      fnm <- filepath %>% basename() %>% stringr::str_remove_all("\\.(laz|las)$")
      new_fnm <- paste0(normalizePath(outdir),"/",fnm,"_epsg",new_crs,".las")
      # reproject
      las <- sf::st_transform(las, paste0("EPSG:", new_crs))
      # write
      lidR::writeLAS(las, file = new_fnm)
      return(new_fnm)
    }
  }
