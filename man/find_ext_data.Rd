% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_ext_data.R
\name{find_ext_data}
\alias{find_ext_data}
\title{Find the location of external data}
\usage{
find_ext_data(
  input_treemap_dir = NULL,
  input_foresttype_dir = NULL,
  input_landfire_dir = NULL
)
}
\arguments{
\item{input_treemap_dir}{character. directory where Treemap 2016 exists. Use \code{\link[=get_treemap]{get_treemap()}} first.}

\item{input_foresttype_dir}{character. directory where Forest Type Groups data exists. Use \code{\link[=get_foresttype]{get_foresttype()}} first.}

\item{input_landfire_dir}{character. directory where LANDFIRE CBD data exists. Use \code{\link[=get_landfire]{get_landfire()}} first.}
}
\value{
Returns a list where the values will be either NULL if unable to locate the external data files
, or the directory where the external data files were located.
The list includes the named variables \code{treemap_dir} and \code{foresttype_dir}
}
\description{
Find the location of external data
Functions \verb{get_*()} download external data
}
\examples{
 \dontrun{
 find_ext_data()
 }

}
