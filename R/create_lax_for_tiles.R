#' @title Create spatial index `.lax` files
#'
#' @param las_file_list a list of .las|.laz files with full directory path
#'
#' @references
#' [https://r-lidar.github.io/lidRbook/spatial-indexing.html](https://r-lidar.github.io/lidRbook/spatial-indexing.html)
#'
#' @description
#' Function to create spatial index files .lax for .las|.laz files to speed up processing
#'
#' @return A list of file names
#'
#' @examples
#'  \dontrun{
#'  f <- list.files(getwd(), pattern = ".*\\.(laz|las)$", full.names = TRUE)
#'  create_lax_for_tiles(las_file_list = f)
#'  }
#' @export
#'
create_lax_for_tiles = function(las_file_list){
    ans <-
      las_file_list %>%
      purrr::map(function(des_file){
        # check validity of data
        l <- rlas::read.lasheader(des_file) %>% length()
        if(l==0 | is.na(l) | is.null(l)){
          return(NULL)
        }else{
          ### Compile the .lax file name
          des_file_lax <- tools::file_path_sans_ext(des_file)
          des_file_lax <- paste0(des_file_lax, ".lax")

          ### See if the .lax version exists in the input directory
          does_file_exist <- file.exists(des_file_lax)

          ### If file does_file_exist, do nothing
          if(does_file_exist == TRUE){return(des_file)}

          ### If file doesnt_file_exist, create a .lax index
          if(does_file_exist == FALSE){
            rlas::writelax(des_file)
            return(des_file)
          }
        }
      }) %>% unlist()
}
