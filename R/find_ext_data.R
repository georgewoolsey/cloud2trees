#' @title Find the location of external data
#'
#' @param input_treemap_dir character. directory where Treemap 2016 exists. Use [get_treemap()] first.
#' @param input_foresttype_dir character. directory where Forest Type Groups data exists. Use [get_foresttype()] first.
#' @param input_landfire_dir character. directory where LANDFIRE CBD data exists. Use [get_landfire()] first.
#'
#' @description
#' Find the location of external data
#' Functions `get_*()` download external data
#'
#' @return Returns a list where the values will be either NULL if unable to locate the external data files
#' , or the directory where the external data files were located.
#' The list includes the named variables `treemap_dir` and `foresttype_dir`
#'
#' @examples
#'  \dontrun{
#'  find_ext_data()
#'  }
#'
#' @export
#'
find_ext_data <- function(
  input_treemap_dir = NULL
  , input_foresttype_dir = NULL
  , input_landfire_dir = NULL
){
  ####################################################
  # let's list all possible directories to check
  ####################################################
    # the package directory from get_url_data()
    pkg_dir <- pkg_dir()

    # user input
    i_treemap_dir <- ifelse(
        !dir.exists( file.path(input_treemap_dir) ) || purrr::is_empty( normalizePath(file.path(input_treemap_dir)) )
        , pkg_dir
        , normalizePath(file.path(input_treemap_dir))
      )
    i_foresttype_dir <- ifelse(
        !dir.exists( file.path(input_foresttype_dir) ) || purrr::is_empty( normalizePath(file.path(input_foresttype_dir)) )
        , pkg_dir
        , normalizePath(file.path(input_foresttype_dir))
      )
    i_landfire_dir <- ifelse(
        !dir.exists( file.path(input_landfire_dir) ) || purrr::is_empty( normalizePath(file.path(input_landfire_dir)) )
        , pkg_dir
        , normalizePath(file.path(input_landfire_dir))
      )

    # current working directory
    pwd_dir <- getwd()

    # stored directory
      # this file is written upon the first successful run of get_treemap
      fp_treemap <- file.path(pkg_dir, "location_treemap.csv")
      if(file.exists(fp_treemap)){
        s_treemap_dir <- fp_treemap %>%
          readr::read_csv(show_col_types = F) %>%
          dplyr::pull(1) %>%
          .[1] %>%
          file.path()
      }else{s_treemap_dir <- pkg_dir}
      # this file is written upon the first successful run of get_foresttype
      fp_foresttype <- file.path(pkg_dir, "location_foresttype.csv")
      if(file.exists(fp_foresttype)){
        s_foresttype_dir <- fp_foresttype %>%
          readr::read_csv(show_col_types = F) %>%
          dplyr::pull(1) %>%
          .[1] %>%
          file.path()
      }else{s_foresttype_dir <- pkg_dir}
      # this file is written upon the first successful run of get_landfire
      fp_landfire <- file.path(pkg_dir, "location_landfire.csv")
      if(file.exists(fp_landfire)){
        s_landfire_dir <- fp_landfire %>%
          readr::read_csv(show_col_types = F) %>%
          dplyr::pull(1) %>%
          .[1] %>%
          file.path()
      }else{s_landfire_dir <- pkg_dir}

  ####################################################
  # data frame the unique directories
  ####################################################
    df <- dplyr::tibble(
        dir = c(
          i_treemap_dir, i_foresttype_dir, i_landfire_dir
          , pkg_dir, pwd_dir
          , s_treemap_dir, s_foresttype_dir, s_landfire_dir
        ) %>% normalizePath()
      ) %>%
      dplyr::distinct()

  ####################################################
  # detect files and return
  ####################################################
    # treemap
    ret_treemap_dir <- df$dir %>%
      purrr::map(\(x) check_dir_files(x, "treemap")) %>%
      purrr::detect(function(x){!is.null(x)}) %>%
      .[1]
    # foresttype
    ret_foresttype_dir <- df$dir %>%
      purrr::map(\(x) check_dir_files(x, "foresttype")) %>%
      purrr::detect(function(x){!is.null(x)}) %>%
      .[1]
    # landfire
    ret_landfire_dir <- df$dir %>%
      purrr::map(\(x) check_dir_files(x, "landfire")) %>%
      purrr::detect(function(x){!is.null(x)}) %>%
      .[1]

  # detect files and return
  return(list(
    treemap_dir = ret_treemap_dir
    , foresttype_dir = ret_foresttype_dir
    , landfire_dir = ret_landfire_dir
  ))

}

####################################################
# function to check dir for required files
####################################################
check_dir_files <- function(dir, which_data = "treemap") {
  if(tolower(which_data)=="treemap"){
    # files to look for
    req_file_list <- list(
      treemap2022()
      , treemap2016()
    )
    # folder to append
    which_data <- tolower(which_data)
  }else if(tolower(which_data)=="foresttype"){
    # files to look for
    req_file_list <- list(
      c("foresttype_lookup.csv", "foresttype.tif")
    )
    # folder to append
    which_data <- tolower(which_data)
  }else if(tolower(which_data)=="landfire"){
    # files to look for
    req_file_list <- list(
      c("lc23_cbd_240.tif")
    )
    # folder to append
    which_data <- tolower(which_data)
  }else{
    return(NULL)
  }

  find_file_list_ans <-
    req_file_list %>%
    # map over the function to find files for each list element
    purrr::map(
      \(x)
      find_file_list(dir = dir, flist = x, which_data = which_data)
    ) %>%
    # keep the first non-null element
    purrr::detect(~ !is.null(.x))

  return(find_file_list_ans)
}
# look in folders for a file list
find_file_list <- function(dir, flist, which_data) {
  # if the dir exists
  if(dir.exists(dir)){
    # check the dir
    f <- list.files(dir) %>% tolower()
    ld <- get_list_diff(flist, f)
    if(length(ld)==0){
      return(dir)
    }else{ # wasn't found in dir
      # check the appended dir
      ad <- file.path(dir, which_data)
      if(dir.exists(ad)){
        f <- list.files(ad) %>% tolower()
        ld <- get_list_diff(flist, f)
      }else{return(NULL)}
      # check if found in appended dir
      if(length(ld)==0){
        return(ad)
      }else{
        return(NULL)
      } # wasn't found in appended dir
    } # wasn't found in dir
  }else{
    return(NULL)
  }
}

####################################################
# function to check dir for tree list data
####################################################
treemap_data_finder <- function(dir) {
  treemap2022 <- treemap2022()
  treemap2016 <- treemap2016()
  f <- list.files(dir) %>% tolower()
  # find matches
  m2022 <- base::intersect(base::unlist(f), base::unlist(treemap2022))
  m2016 <- base::intersect(base::unlist(f), base::unlist(treemap2016))
  if(length(m2022)==length(treemap2022)){
    return(list(
      which_treemap = 2022
      , treemap_rast = file.path( dir, treemap2022[1] )
      , treemap_trees = file.path( dir, treemap2022[2] )
    ))
  }else if(length(m2016)==length(treemap2016)){
    return(list(
      which_treemap = 2016
      , treemap_rast = file.path( dir, treemap2016[1] )
      , treemap_trees = file.path( dir, treemap2016[2] )
    ))
  }else{
    stop("couldn't find any treemap data :[")
  }
}
