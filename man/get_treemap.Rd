% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_treemap.R
\name{get_treemap}
\alias{get_treemap}
\title{Download TreeMap 2016 data}
\usage{
get_treemap(savedir = NULL, force = F)
}
\arguments{
\item{savedir}{Optional directory to save data in a new location. Defaults to package contents.}

\item{force}{Whether to overwrite existing data}
}
\description{
To estimate individual tree DBH based on the point cloud detected tree height models are fit to FIA plot data within a buffer of the point cloud boundary.
FIA plots are identified using TreeMap 2016, a model of FIA plot locations imputed throughout forested areas of the conterminous United States at 30 m spatial resolution.
See \code{\link[=trees_dbh]{trees_dbh()}}
}
\examples{
 \dontrun{
 get_treemap()
 }
}
\references{
\url{https://doi.org/10.2737/RDS-2021-0074}
Riley, Karin L.; Grenfell, Isaac C.; Finney, Mark A.; Shaw, John D. 2021. TreeMap 2016: A tree-level model of the forests of the conterminous United States circa 2016. Fort Collins, CO: Forest Service Research Data Archive.
}
