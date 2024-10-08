% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lasr_chm.R
\name{lasr_chm}
\alias{lasr_chm}
\title{Create CHM and apply pits and spikes filling via \code{lasR}}
\usage{
lasr_chm(
  chm_file_name,
  chm_res = 0.25,
  min_height_m = 2,
  max_height_m = 70,
  lap_sz = 3
)
}
\arguments{
\item{chm_file_name}{string. Where to write the CHM.}

\item{chm_res}{numeric. The desired resolution of the CHM produced in meters.}

\item{min_height_m}{numeric. Set the minimum height (m) for individual tree detection}

\item{max_height_m}{numeric. Set the maximum height (m) for the canopy height model}

\item{lap_sz}{numeric. Size of the Laplacian filter kernel (integer value, in pixels) for \code{\link[lasR:pit_fill]{lasR::pit_fill()}}}
}
\value{
A \code{lasR} pipeline
}
\description{
use \code{lasR} to create CHM and apply pits and spikes filling for raster based on St-Onge 2008 (see reference).
}
\references{
\url{https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=81365288221f3ac34b51a82e2cfed8d58defb10e}
}
\keyword{internal}
