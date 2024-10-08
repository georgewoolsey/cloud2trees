% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_lax_for_tiles.R
\name{create_lax_for_tiles}
\alias{create_lax_for_tiles}
\title{Create spatial index \code{.lax} files}
\usage{
create_lax_for_tiles(las_file_list)
}
\arguments{
\item{las_file_list}{a list of .las|.laz files with full directory path}
}
\value{
A list of file names
}
\description{
Function to create spatial index files .lax for .las|.laz files to speed up processing
}
\examples{
 \dontrun{
 f <- list.files(getwd(), pattern = ".*\\\\.(laz|las)$", full.names = TRUE)
 create_lax_for_tiles(las_file_list = f)
 }
}
\references{
\url{https://r-lidar.github.io/lidRbook/spatial-indexing.html}
}
