# Use a CHM raster to detect individual trees

`raster2trees()` is an all-in-one function to process a CHM raster and
return a spatial data frame of tree crown polygons and points. The order
of operations is:

- Perform individual tree detection using
  [`lidR::locate_trees()`](https://rdrr.io/pkg/lidR/man/locate_trees.html)
  with the [`lidR::lmf()`](https://rdrr.io/pkg/lidR/man/itd_lmf.html)
  algorithm

- Delineate tree crowns using
  [`ForestTools::mcws()`](https://rdrr.io/pkg/ForestTools/man/mcws.html)

Note, this function does not estimate DBH for the detected trees and
only returns tree location, crown area, and height information. To
estimate tree DBH from the detected tree heights see
[`trees_dbh()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_dbh.md).

## Usage

``` r
raster2trees(
  chm_rast,
  outfolder,
  ws = itd_ws_functions()[["log_fn"]],
  min_height = 2,
  min_crown_area = 0.1,
  tempdir = tempdir()
)
```

## Arguments

- chm_rast:

  raster. A raster from `terra` or `stars`representing a canopy height
  model

- outfolder:

  string. The path of a folder to write the crown vector data to

- ws:

  numeric or function. Length or diameter of the moving window used to
  detect the local maxima in the units of the input data (usually
  meters). If it is numeric a fixed window size is used. If it is a
  function, the function determines the size of the window at any given
  location on the canopy. By default function takes the height of a
  given pixel as its only argument and return the desired size of the
  search window when centered on that pixel.

- min_height:

  numeric. Set the minimum height (m) for individual tree detection

- min_crown_area:

  numeric. Set the minimum crown area (m2) for individual tree detection

- tempdir:

  string. Directory to write intermediate files. Intermediate files are
  only created for large rasters too big to fit in memory.

## Value

Returns a spatial data frame of individual tree crown vectors detected
using the CHM. The tree top point coordinates are located in the
`tree_x` and `tree_y` columns. The process also writes two `.gpkg` files
to the `outfolder` directory: `chm_detected_crowns.gpkg` and
`chm_detected_tree_tops.gpkg`

## References

<https://r-lidar.github.io/lidRbook/itd.html>

## Examples

``` r
 if (FALSE) { # \dontrun{
 f <- paste0(system.file(package = "cloud2trees"),"/extdata/chm.tif")
 crowns_sf <- raster2trees(chm_rast = terra::rast(f), outfolder = tempdir())
 crowns_sf %>% class()
 crowns_sf %>% dplyr::glimpse()
 crowns_sf %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(fill=tree_height_m))
 } # }
```
