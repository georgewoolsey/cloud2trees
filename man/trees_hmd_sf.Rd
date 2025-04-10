% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_hmd.R
\name{trees_hmd_sf}
\alias{trees_hmd_sf}
\title{Estimate HMD using tree crown polygons and normalized point cloud data}
\usage{
trees_hmd_sf(
  trees_poly,
  norm_las = NULL,
  tree_sample_n = NA,
  tree_sample_prop = NA,
  force_same_crs = F,
  trees_sample = NA,
  ofile = NA
)
}
\arguments{
\item{trees_poly}{sf. A \code{sf} class object with POLYGON geometry (see \code{\link[sf:st_geometry_type]{sf::st_geometry_type()}}), the program will use the data "as-is" and only require the \code{treeID} and \code{tree_height_m} columns.
Or the path to a single spatial polygon file.}

\item{norm_las}{character. a directory with nomalized las files, the path of a single .laz|.las file", -or- an object of class \code{LAScatalog}.
It is your responsibility to ensure that the point cloud is projected the same as the \code{trees_poly} data}

\item{tree_sample_n, tree_sample_prop}{numeric. Provide either \code{tree_sample_n}, the number of trees, or \code{tree_sample_prop}, the
proportion of the trees to attempt to extract a HMD from the point cloud for.
If neither are supplied, \code{tree_sample_n = 777} will be used. If both are supplied, \code{tree_sample_n} will be used.
Increasing \code{tree_sample_prop} toward one (1) will increase the processing time, perhaps significantly depending on the number of trees in the \code{trees_poly} data.}

\item{force_same_crs}{logical. force the same crs between the point cloud and polygon if confident that data are in same projection.
data created by a \code{cloud2trees} pipeline (e.g. \code{\link[=cloud2raster]{cloud2raster()}}) will always have the same projection even if not recognized by \code{lidR} functions}

\item{trees_sample}{data.frame. provide your own tree sample list such as one generated from \code{sample_trees_flist()} that includes the \code{treeID} column.
If provided, the tree_sample_n,tree_sample_prop will be ignored}

\item{ofile}{character or logical. if a character value is provided the output will be written to the disk as a csv at the location provided.
If set to TRUE and a file path was used as the input for \code{trees_poly}, then a csv file will be written to the same location with the same name prefixed with "hmd_".
Leave as NA to return a data.frame of the trees from tree list from \code{trees_poly} with HMD values added}
}
\description{
\code{trees_hmd_sf()} uses the input tree crown polygons (e.g. as exported by \code{\link[=raster2trees]{raster2trees()}}) with the columns
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
\keyword{internal}
