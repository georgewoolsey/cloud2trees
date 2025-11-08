# Individual Tree Detection (ITD) tuning

`itd_tuning()` is used to visually assess tree crown delineation results
from different window size functions used for the detection of
individual trees. The `cloud2trees` package performs individual tree
detection using
[`lidR::locate_trees()`](https://rdrr.io/pkg/lidR/man/locate_trees.html)
with the [`lidR::lmf()`](https://rdrr.io/pkg/lidR/man/itd_lmf.html)
algorithm. The local maximum filter algorithm allows for a constant
window size or a variable window size defined by a function. See the
`lidR` [package book](https://r-lidar.github.io/lidRbook/itd.html) for
excellent detail on ITD and defining window size.

`itd_tuning()` allows users to test different window size functions on a
sample of data to determine which function is most suitable for the area
being analyzed. The preferred function can then be used in the `ws`
parameter in
[`raster2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/raster2trees.md)
and/or
[`cloud2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/cloud2trees.md).

## Usage

``` r
itd_tuning(
  input_las_dir,
  n_samples = 3,
  ws_fn_list = NULL,
  min_height = 2,
  chm_res_m = 0.25
)
```

## Arguments

- input_las_dir:

  character. directory where .las\|.laz point cloud data
  exists...program will search all sub-directories for all .las\|.laz
  files and process them as one

- n_samples:

  numeric. The number of sample plots of 0.1 ha on which to test the
  window functions. The maximum is 5. The center of the point cloud data
  coverage will always be the first plot sampled so long as points exist
  in the central 0.1 ha.

- ws_fn_list:

  list. A function or a named list of functions. Leave as NULL to test
  default exponential (concave up), linear, and logarithmic (concave
  down) functions. If providing a custom function, it must always return
  a numeric value \>0 (see examples).

- min_height:

  numeric. Set the minimum height (m) for individual tree detection

- chm_res_m:

  numeric. The desired resolution of the CHM produced in meters.

## Value

Returns a list with: 1) "plot_samples" is a plot of the sample canopy
height model (CHM) and extracted tree crowns for each window size
tested; and 2) "ws_fn_list" is a list of the window size functions
tested which can be used to pass the desired function to the `ws`
parameter in
[`raster2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/raster2trees.md)
and/or
[`cloud2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/cloud2trees.md).

## References

<https://r-lidar.github.io/lidRbook/itd.html>

## Examples

``` r
 if (FALSE) { # \dontrun{
  # do it
  library(tidyverse)
  # test las file but this could also be a directory path with >1 .las|.laz files
  i <- system.file(package = "lidR", "extdata", "MixedConifer.laz")
  ####################################################
  # check the default itd_tuning() window functions
  ####################################################
   # run it with defaults
   itd_tuning_ans <- itd_tuning(input_las_dir = i)
   # what's in it?
   names(itd_tuning_ans)
   # look at the tuning plot showing the tree crowns on the CHM
   itd_tuning_ans$plot_samples
   # look at the summary of the trees detected by each ITD function
   itd_tuning_ans$plot_sample_summary
   # the "exp_fn" looks pretty good, let's store it
   best_default <- itd_tuning_ans$ws_fn_list$exp_fn
   # we can see what this function looks like for window size
   ggplot2::ggplot() +
     ggplot2::geom_function(fun = best_default) +
     ggplot2::xlim(-5,60) +
     ggplot2::labs(x = "heights", y = "ws", color = "")
   # pass our best function to the cloud2trees() to process the full point cloud coverage
   cloud2trees_ans <- cloud2trees(output_dir = tempdir(), input_las_dir = i, ws = best_default)
   # the same plot as the the tuning plot with tree crowns overlaid on CHM
   ggplot2::ggplot() +
     ggplot2::geom_tile(
       data = cloud2trees_ans$chm_rast %>%
         terra::as.data.frame(xy=T) %>%
         dplyr::rename(f=3)
       , mapping = ggplot2::aes(x = x, y = y, fill = f)
       , na.rm = T
     ) +
     ggplot2::scale_fill_viridis_c(
       option = "plasma"
       , breaks = scales::breaks_extended(n=10)
     ) +
     ggplot2::geom_sf(
       data = cloud2trees_ans$crowns_sf
       , fill = NA, color = "gray33", lwd = 1
     ) +
     ggplot2::scale_x_continuous(expand = c(0, 0)) +
     ggplot2::scale_y_continuous(expand = c(0, 0)) +
     ggplot2::labs(x = "", y = "", fill = "CHM (m)") +
     ggplot2::theme_light() +
     ggplot2::theme(axis.text = ggplot2::element_blank())
  ####################################################
  # let's test some custom window functions
  ####################################################
    # a constant window size has to be defined as:
     ## x*0 + constant
     my_constant <- function(x){(x * 0) + 3} ## will always return 3
    # a custom linear function
     my_linear <- function(x) {(x * 0.1) + 3}
    # run it with custom functions
     itd_tuning_ans2 <- itd_tuning(
       input_las_dir = i
       , ws_fn_list = list(
          my_constant=my_constant
          , my_linear=my_linear
          , best_default=best_default # the best from our first test
        )
       , n_samples = 2
      )
    # look at the tuning plot showing the tree crowns on the CHM
     itd_tuning_ans$plot_samples
    # look at the summary of the trees detected by each ITD function
     itd_tuning_ans$plot_sample_summary
    # we can see what our custom "my_linear" function looks like
     ggplot2::ggplot() +
       ggplot2::geom_function(fun = itd_tuning_ans2$ws_fn_list$my_linear) +
       ggplot2::xlim(-5,60) +
       ggplot2::labs(x = "heights", y = "ws", color = "")
 } # }
```
