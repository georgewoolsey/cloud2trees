% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/itd_ws_functions.R
\name{itd_ws_functions}
\alias{itd_ws_functions}
\title{Individual Tree Detection (ITD) functions}
\usage{
itd_ws_functions()
}
\value{
Returns a list of named functions which can be used to pass the desired function to the \code{ws}
parameter in \code{\link[=raster2trees]{raster2trees()}} and/or \code{\link[=cloud2trees]{cloud2trees()}}.
}
\description{
\code{itd_ws_functions()} is a list of functions that can be used for determining a
variable window size for the detection of individual trees.
The \code{cloud2trees} package performs individual tree detection using \code{\link[lidR:locate_trees]{lidR::locate_trees()}} with the \code{\link[lidR:itd_lmf]{lidR::lmf()}} algorithm.
The local maximum filter algorithm allows for a constant window size or a variable window size defined by a function.
See the \code{lidR} \href{https://r-lidar.github.io/lidRbook/itd.html}{package book} for excellent detail on ITD and defining window size.
}
\examples{
 \dontrun{
  # what is this?
  itd_ws_functions() \%>\% class()
  # is this list named?
  itd_ws_functions() \%>\% names()
  # what is the first thing in the list named?
  itd_ws_functions()[1] \%>\% names()
  # we can reference it by name
  itd_ws_functions()["lin_fn"] \%>\% names()
  # how can we access a function?
  itd_ws_functions()["exp_fn"] \%>\% is.function() # still a list
  itd_ws_functions()[["exp_fn"]] \%>\% is.function() # now a function
  itd_ws_functions() \%>\% purrr::pluck("exp_fn") \%>\% is.function() # also now a function
  # let's store it
  def_exp_fn <- itd_ws_functions()[["exp_fn"]]
  # can we use it?
  def_exp_fn(9)
  # can we plot a function?
  ggplot2::ggplot() +
    ggplot2::geom_function(fun = itd_ws_functions()[["exp_fn"]]) +
    ggplot2::xlim(-1,60)
  # can we plot all functions?
  ggplot2::ggplot() +
    ggplot2::geom_function(ggplot2::aes(color="lin_fn"), fun = itd_ws_functions()[["lin_fn"]]) +
    ggplot2::geom_function(ggplot2::aes(color="exp_fn"), fun = itd_ws_functions()[["exp_fn"]]) +
    ggplot2::geom_function(ggplot2::aes(color="log_fn"), fun = itd_ws_functions()[["log_fn"]]) +
    ggplot2::xlim(-1,60)
 }
}
