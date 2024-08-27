#' @title Download TreeMap 2016 data
#' @param force Whether to overwrite exising data
#' @param savedir Optional directory to save data in a new location. Defaults to package contents.
#' @references
#' Riley, Karin L.; Grenfell, Isaac C.; Finney, Mark A.; Shaw, John D. 2021. TreeMap 2016: A tree-level model of the forests of the conterminous United States circa 2016. Fort Collins, CO: Forest Service Research Data Archive. https://doi.org/10.2737/RDS-2021-0074
#'
#' @description
#' To estimate individual tree DBH based on the point cloud detected tree height models are fit to FIA plot data within a buffer of the point cloud boundary.
#' FIA plots are identified using TreeMap 2016, a model of FIA plot locations imputed throughout forested areas of the conterminous United States at 30 m spatial resolution.
#' @examples
#'  \donttest{
#'  get_treemap()
#'  }
#' @export
#'
get_treemap <- function(savedir=NULL,force=F){
  #Store users timeout options
  timeout_option_backup <- getOption("timeout")
  options(timeout = max(3600, getOption("timeout")))

  if(is.null(savedir)){
    # create dir
    dir.create(paste0(system.file(package = "cloud2trees"),"/extdata"), showWarnings = FALSE)
    # names
    destination <- paste0(system.file(package = "cloud2trees"),"/extdata/treemap.zip")
    dirname <- paste0(system.file(package = "cloud2trees"),"/extdata/treemap")
  } else{
    destination <- file.path(savedir,"treemap.zip")
    dirname <- file.path(savedir,"treemap")
  }
  # create dir
  dir.create(dirname, showWarnings = FALSE)

  #check if already exists.
  f_dir <- file.path(dirname(destination),"treemap/")
  # f_dir <- paste0(system.file("extdata", "treemap/", package = "cloud2trees"))
  f <- toupper(list.files(f_dir))
  if(length(f)==0){f <- ""}
  if(
    max(grepl("TREEMAP2016.TIF", f))==1 & max(grepl("TREEMAP2016_TREE_TABLE.CSV", f))==1
  ){
    if(!force){
      warning(paste("Data has already been downloaded to",dirname,", use force=T to overwrite"))
      return(NULL)
    }
  }

  # get data
  eval_url <- "https://s3-us-west-2.amazonaws.com/fs.usda.rds/RDS-2021-0074/RDS-2021-0074_Data.zip"
  message(paste("Downloading file to",destination))
  download.file(eval_url, destination, mode = "wb")
  unzip_download(destination)

  options(timeout = timeout_option_backup)
}
## unzip function
unzip_download <- function(destination){
  #location of unzip
  base_dir <- dirname(destination)

  #get file names
  unzip_folder <- unzip(destination, list = TRUE)$Name[1]
  unzipped_folder <- file.path(base_dir,unzip_folder)
  unzip(destination,exdir=base_dir)
  final_name <- file.path(base_dir,"treemap/")

  #Force delete of any previous folder
  unlink(final_name,recursive = T)
  file.rename(unzipped_folder,final_name)

  #Remove zipped files
  unlink(destination)
}
