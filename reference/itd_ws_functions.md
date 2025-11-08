# Individual Tree Detection (ITD) functions

`itd_ws_functions()` is a list of functions that can be used for
determining a variable window size for the detection of individual
trees. The `cloud2trees` package performs individual tree detection
using
[`lidR::locate_trees()`](https://rdrr.io/pkg/lidR/man/locate_trees.html)
with the [`lidR::lmf()`](https://rdrr.io/pkg/lidR/man/itd_lmf.html)
algorithm. The local maximum filter algorithm allows for a constant
window size or a variable window size defined by a function. See the
`lidR` [package book](https://r-lidar.github.io/lidRbook/itd.html) for
excellent detail on ITD and defining window size.

## Usage

``` r
itd_ws_functions()
```

## Value

Returns a list of named functions which can be used to pass the desired
function to the `ws` parameter in
[`raster2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/raster2trees.md)
and/or
[`cloud2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/cloud2trees.md).

## Examples

``` r
 if (FALSE) { # \dontrun{
  # what is this?
  itd_ws_functions() %>% class()
  # is this list named?
  itd_ws_functions() %>% names()
  # what is the first thing in the list named?
  itd_ws_functions()[1] %>% names()
  # we can reference it by name
  itd_ws_functions()["lin_fn"] %>% names()
  # how can we access a function?
  itd_ws_functions()["exp_fn"] %>% is.function() # still a list
  itd_ws_functions()[["exp_fn"]] %>% is.function() # now a function
  itd_ws_functions() %>% purrr::pluck("exp_fn") %>% is.function() # also now a function
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
 } # }
```
