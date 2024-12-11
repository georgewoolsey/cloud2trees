#' @title Download Forest Type Groups of the Continental United States data
#' @param force Whether to overwrite exising data
#' @param savedir Optional directory to save data in a new location. Defaults to package contents.
#' @references
#' * [Forest Type Groups of the Continental United States](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d)
#' Wilson, B.T. (2023). Forest Type Groups of the Continental United States.
#'
#' @description
#' The Forest Type Groups of the Continental United States data is used to estimate individual tree forest type group.
#' See [trees_species()]
#'
#' @examples
#'  \dontrun{
#'  get_foresttype()
#'  }
#' @export
#'
get_foresttype <- function(
  savedir = NULL
  , force = F
){
  # set up parameters to pass to get_url_data()
  my_eval_url <- "https://zenodo.org/records/14343811/files/foresttype.zip?download=1"
  my_my_name <- "foresttype"
  my_req_file_list <- c("foresttype_lookup.csv", "foresttype.tif")
  my_cleanup_zip <- T

  # call get_url_data()
  get_url_data(
    eval_url = my_eval_url
    , my_name = my_my_name
    , savedir = savedir
    , req_file_list = my_req_file_list
    , force = force
    , cleanup_zip = my_cleanup_zip
  )
}
