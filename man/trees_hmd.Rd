% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trees_hmd.R
\name{trees_hmd}
\alias{trees_hmd}
\title{Estimate HMD using tree crown polygons and normalized point cloud data}
\usage{
trees_hmd(
  trees_poly,
  norm_las = NULL,
  tree_sample_n = NA,
  tree_sample_prop = NA,
  estimate_missing_hmd = TRUE,
  force_same_crs = F,
  outfolder = tempdir()
)
}
\arguments{
\item{trees_poly}{must be one of the following that has required attributes \code{treeID} and \code{tree_height_m}:
\itemize{
\item \code{sf} class object with POLYGON geometry (see \code{\link[sf:st_geometry_type]{sf::st_geometry_type()}}). Recommended for smaller tree lists (e.g. <100k) that can fit in memory.
\item character vector with the path to a single or multiple spatial files that can be read by \code{\link[sf:st_read]{sf::st_read()}} and have POLYGON geometry. Recommended for large tree lists (e.g. 100k+) that might cause memory issues.
\item character with the path to a directory that has "final_detected_crowns*" files from \code{\link[=cloud2trees]{cloud2trees()}} or \code{\link[=raster2trees]{raster2trees()}}. Recommended for large tree lists (e.g. 100k+) that might cause memory issues.
}}

\item{norm_las}{character. a directory with nomalized las files, the path of a single .laz|.las file", -or- an object of class \code{LAScatalog}.
It is your responsibility to ensure that the point cloud is projected the same as the \code{trees_poly} data}

\item{tree_sample_n, tree_sample_prop}{numeric. Provide either \code{tree_sample_n}, the number of trees, or \code{tree_sample_prop}, the
proportion of the trees to attempt to extract a HMD from the point cloud for.
If neither are supplied, \code{tree_sample_n = 777} will be used. If both are supplied, \code{tree_sample_n} will be used.
Increasing \code{tree_sample_prop} toward one (1) will increase the processing time, perhaps significantly depending on the number of trees in the \code{trees_poly} data.}

\item{estimate_missing_hmd}{logical. it is not likely that HMD will be extracted successfully from every tree (especially in low density clouds).
Should the missing HMD values be estimated using the tree height and location information based on trees for which HMD is successfully extracted?}

\item{force_same_crs}{logical. force the same crs between the point cloud and polygon if confident that data are in same projection.
data created by a \code{cloud2trees} pipeline (e.g. \code{\link[=cloud2raster]{cloud2raster()}}) will always have the same projection even if not recognized by \code{lidR} functions}

\item{outfolder}{string. The path of a folder to write the model data to.
Note, in the actual missing value estimation many RF models are estimated and model averaging is used.
However, only the first estimated model is saved in this export which does not fully represent the process used to fill in missing values.}
}
\value{
Returns a spatial data frame of individual trees with the added columns: \code{max_crown_diam_height_m}, \code{is_training_hmd}
}
\description{
\code{trees_hmd()} uses the input tree crown polygons (e.g. as exported by \code{\link[=raster2trees]{raster2trees()}}) with the columns
\code{treeID} and \code{tree_height_m} to extracting the height of the maximum crown diameter (HMD) using height normalized point cloud data (e.g. as exported by \code{\link[=cloud2raster]{cloud2raster()}}).

HMD is extracted directly from the height normalized point cloud by finding the height of the non-ground point farthest from the tree center (i.e. tree top).

An early version of this process was developed by \href{https://github.com/bi0m3trics}{Andrew Sanchez Meador}.

There are likely to be trees for which there is insufficient data in the point cloud to successfully estimate HMD. The user can elect to estimate missing HMD values which is accomplished via:
\itemize{
\item Attempt to extract HMD from all trees
\item Successfully extracted HMD trees become training data used to estimate the height-HMD allometry relationship that is spatially informed using the relative tree location compared to the training data
\item The height and location predicting HMD model built from the point cloud training data is used to predict HMD for the non-training (i.e. missing HMD) data
}
}
\examples{
 \dontrun{
  library(tidyverse)
  library(sf)
  # example tree crown polygons
  f <- paste0(system.file(package = "cloud2trees"),"/extdata/crowns_poly.gpkg")
  crowns <- sf::st_read(f, quiet = T)
  # example normalized las files are in this directory
  norm_d <- paste0(system.file(package = "cloud2trees"),"/extdata/norm_las")
  # now run the trees_hmd()
  trees_hmd_ans <- trees_hmd(
     trees_poly = crowns
     , norm_las = norm_d
     , force_same_crs = T
     , estimate_missing_hmd = T)
  # what?
  trees_hmd_ans \%>\% dplyr::glimpse()
  # spatial polygons
  trees_hmd_ans \%>\% ggplot2::ggplot() +
     ggplot2::geom_sf(ggplot2::aes(fill=max_crown_diam_height_m))
  # relationship between height and hmd
  trees_hmd_ans \%>\%
     ggplot2::ggplot(
       ggplot2::aes(
         x = tree_height_m, y = max_crown_diam_height_m, color=is_training_hmd
       )
     ) +
     ggplot2::geom_point()
  #### try a file list
  #### Recommended for large tree lists (e.g. 100k+) that might cause memory issues.
  # we'll split the crowns
  # as is done automatically for tree lists >250k by raster2trees() and cloud2trees()
  crowns <- crowns \%>\%
    dplyr::mutate(
      # makes 2 groups of data
      grp = ceiling(dplyr::row_number()/(dplyr::n()/2))
    )
  # make file names
  my_dir <- tempdir()
  fnm_1 <- file.path(my_dir, "crowns1.gpkg")
  fnm_2 <- file.path(my_dir, "crowns2.gpkg")
  fnm_1
  # write the data
  sf::st_write(crowns \%>\% dplyr::filter(grp==1), dsn = fnm_1, append = F) # grp 1
  sf::st_write(crowns \%>\% dplyr::filter(grp==2), dsn = fnm_2, append = F) # grp 2
  # try trees_cbh with our file list
  flist <- c(fnm_1,fnm_2)
  # now run the trees_hmd()
  trees_hmd_ans2 <- trees_hmd(
     trees_poly = flist
     , norm_las = norm_d
     , force_same_crs = T
     , estimate_missing_hmd = T)
  # tabulate training data
  trees_hmd_ans \%>\%
    sf::st_drop_geometry() \%>\%
    dplyr::count(is_training_hmd)
 }
}
\references{
An early version of this process was developed by \href{https://github.com/bi0m3trics}{Andrew Sanchez Meador}.
}
