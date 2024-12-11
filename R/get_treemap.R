#' @title Download TreeMap 2016 data
#' @param force Whether to overwrite exising data
#' @param savedir Optional directory to save data in a new location. Defaults to package contents.
#' @references
#' [https://doi.org/10.2737/RDS-2021-0074](https://doi.org/10.2737/RDS-2021-0074)
#' Riley, Karin L.; Grenfell, Isaac C.; Finney, Mark A.; Shaw, John D. 2021. TreeMap 2016: A tree-level model of the forests of the conterminous United States circa 2016. Fort Collins, CO: Forest Service Research Data Archive.
#'
#' @description
#' To estimate individual tree DBH based on the point cloud detected tree height models are fit to FIA plot data within a buffer of the point cloud boundary.
#' FIA plots are identified using TreeMap 2016, a model of FIA plot locations imputed throughout forested areas of the conterminous United States at 30 m spatial resolution.
#' See [trees_dbh()]
#'
#' @examples
#'  \dontrun{
#'  get_treemap()
#'  }
#' @export
#'
get_treemap <- function(
  savedir = NULL
  , force = F
){
  # set up parameters to pass to get_url_data()
    # old url: "https://s3-us-west-2.amazonaws.com/fs.usda.rds/RDS-2021-0074/RDS-2021-0074_Data.zip"
  my_eval_url <- "https://usfs-public.box.com/shared/static/yz7h8b8v92scoqfwukjyulokaevzo6v6.zip" # updated 2024-12-10
  my_my_name <- "treemap"
  my_req_file_list <- c("treemap2016.tif", "treemap2016_tree_table.csv")
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
