#' @title Download url data
#' @param eval_url Required url of the .zip file to be downloaded. url string must end in .zip
#' @param my_name Required name for the folder which data will be extracted to
#' @param savedir Optional directory to save data in a new location. Defaults to package contents.
#' @param req_file_list Optional list of files to check for before re-downloading full data
#' @param force Whether to overwrite exising data
#' @param cleanup_zip Whether to remove the .zip file after extracting the contents
#'
#' @description
#' Generic function to download data from a url with .zip file data
#' Functions [get_treemap()] and [get_foresttype()] use this
#'
#' @examples
#'  \dontrun{
#'  get_treemap()
#'  }
#'
#' @keywords internal
#'
get_url_data <- function(
  eval_url = NULL
  , my_name = NULL
  , savedir = NULL
  , req_file_list = NULL
  , force = F
  , cleanup_zip = T
  , move_files_to_top = T
){
  if(is.null(my_name)){stop("must provide a name for the folder which data will be extracted to in the my_name parameter")}
  # check url
  has_zip <- stringr::word(dplyr::coalesce(eval_url, "heyxxxxx"), -1, sep = "/") %>%
    tolower() %>%
    stringr::str_detect(".zip")
  if(is.null(eval_url) || !inherits(eval_url, "character") || has_zip==F){
    stop("must provide a url character string in the eval_url parameter that results in the download of a .zip file")
  }

  #Store users timeout options
  timeout_option_backup <- getOption("timeout")
  options(timeout = max(3600, getOption("timeout")))

  ## directory names
  if(is.null(savedir)){
    # create dir
    dir.create(file.path(system.file(package = "cloud2trees"),"extdata"), showWarnings = FALSE)
    # names
    dirname <- file.path(system.file(package = "cloud2trees"),"extdata", my_name)
  }else{
    dirname <- file.path(savedir,my_name)
  }
  # create dir
  dir.create(dirname, showWarnings = FALSE)
  # full file zip path
  zip_file_pth_nm <- file.path(dirname, paste0(my_name,".zip"))

  # check if required files already exists.
  if(is.null(req_file_list)){req_file_list <- "hey_xxxxxx"}
  f <- list.files(dirname) %>% tolower()
  ld <- get_list_diff(req_file_list, f)
  if(length(ld)==0){
    if(!force){
      warning(paste("Data has already been downloaded to",dirname,", use force=T to overwrite"))
      return(NULL)
    }
  }

  # get url data
  message(paste("Downloading file to",zip_file_pth_nm))
  download.file(eval_url, zip_file_pth_nm, mode = "wb")

  # unzip
  unzip_download(zip_file_pth_nm, zip_rm = cleanup_zip, move_to_top = move_files_to_top)

  options(timeout = timeout_option_backup)
}
####################################################################
## intermediate fn to get list difference
####################################################################
get_list_diff <- function(x, y) {
  if(inherits(x, "character") & inherits(y, "character")){
      d <- x[!(x %in% y)]
      d <- unique(d)
  }else{
    stop("must provide character list in x and y")
  }
  return(d)
}
####################################################################
## intermediate fn to unzip
####################################################################
unzip_download <- function(destination, zip_rm = TRUE, move_to_top = TRUE){
  #location of unzip
  base_dir <- dirname(destination)

  #get file and directory names
  uz <- unzip(destination, exdir=base_dir, list = T)

  # keep only list of files in a child directory
  uz <- uz %>%
      dplyr::filter(
        !stringr::str_ends(Name, "/") # is not a directory
        & stringr::str_detect(Name, "/") # is in child directory
      )

  # unzip it
  unzip(destination, exdir=base_dir)

  # if requested to move files to top folder
  if(move_to_top==T && dplyr::coalesce(nrow(uz),0)>0){
    # move all files
    uz$Name %>%
    purrr::map(function(x){
      # from file full path
      f <- file.path(base_dir, x)
      # to file full path
      t <- file.path(base_dir, basename(x))
      # copy
      file.copy(f, to = t, overwrite = T)
      # remove original
      file.remove(f)
    })
    # # remove empty directories
    # uz$Name %>%
    #   dirname() %>%
    #   unique() %>%
    #   purrr::map(function(x){
    #     # full directory
    #     d <- file.path(base_dir, x)
    #     # remove original
    #     file.remove(d)
    #   })
  }

  #delete the zip file
  if(zip_rm==T){
    unlink(destination)
  }
}
