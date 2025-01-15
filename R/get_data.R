#' @title Download all external data used by package
#' @param savedir Optional directory to save data in a new location. Defaults to package contents.
#' @param force Whether to overwrite existing data
#'
#' @description
#' The package requires external data to estimate individual tree DBH (see [get_treemap()]) and to extract the forest type for a tree list (see [get_foresttype()]).
#' This is a all-in-one function that downloads all of the external data used by the package.
#'
#' @examples
#'  \dontrun{
#'  get_data()
#'  }
#' @export
#'
get_data <- function(
  savedir = NULL
  , force = F
){
  # call individual download functions
  get_treemap(savedir = savedir, force = force)
  get_foresttype(savedir = savedir, force = force)
  get_landfire(savedir = savedir, force = force)
}
