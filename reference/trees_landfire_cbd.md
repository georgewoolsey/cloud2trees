# Extract LANDFIRE CBD raster cell value for a tree list based on location

`trees_landfire_cbd()` uses the input tree list (e.g. as exported by
[`raster2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/raster2trees.md))
with the columns `treeID`, `tree_x`, `tree_y` to attach LANDFIRE's
Forest Canopy Bulk Density (CBD) data estimate in kilograms per cubic
meter produced jointly by the U.S. Department of Agriculture and U.S.
Department of the Interior. If a spatial data frame of points is the
input tree list, then the columns `tree_x`, `tree_y` are not required.

LANDFIRE's Forest Canopy Bulk Density (CBD) data is attached to each
tree in the tree list based on the spatial overlap with the raster data
set (see references). Canopy Bulk Density is mass of flammable material
per unit volume of the tree crown typically expressed in units of mass
per unit volume (e.g., kilograms per cubic meter).

The simplified process for attaching the LANDFIRE CBD raster cell value
to a tree is:

- Nearest neighbor imputation is used to fill LANDFIRE data if a tree
  falls inside a non-forest cell in the original data

- The LANDFIRE raster cell value in kilograms per cubic meter is applied
  to a tree based on spatial overlap

## Usage

``` r
trees_landfire_cbd(
  tree_list,
  crs = NA,
  study_boundary = NA,
  input_landfire_dir = NULL,
  max_search_dist_m = 1000
)
```

## Arguments

- tree_list:

  data.frame. A data frame with the columns `treeID`, `tree_x`,
  `tree_y`, and `tree_height_m`. If an `sf` class object with POINT
  geometry (see
  [`sf::st_geometry_type()`](https://r-spatial.github.io/sf/reference/st_geometry_type.html)),
  the program will use the data "as-is" and only require the `treeID`
  column.

- crs:

  string. A crs string as returned from
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)
  or the EPSG code of the x,y coordinates. Defaults to the crs of the
  `tree_list` data if of class "sf".

- study_boundary:

  sf. The boundary of the study are to define the area of the regional
  model. If no boundary given, regional model will be built from
  location of trees in the tree list.

- input_landfire_dir:

  directory where LANDFIRE CBD data exists. Use
  [`get_landfire()`](https://georgewoolsey.github.io/cloud2trees/reference/get_landfire.md)
  first.

- max_search_dist_m:

  number. Maximum search distance (m) to obtain forest type group data
  for trees in `tree_list` that overlap with non-forest data in the
  original Wilson (2023) data. Larger search distances will increase
  processing time and possibly result in memory issues.

## Value

Returns a list of objects: tree_list = spatial data frame of individual
trees with the column `landfire_cell_kg_per_m3` added ; landfire_rast =
raster of kilograms per cubic meter in the area.

## References

- [LANDFIRE Forest Canopy Bulk Density
  (CBD)](https://landfire.gov/fuel/cbd) U.S. Department of Agriculture
  and U.S. Department of the Interior.

## Examples

``` r
 if (FALSE) { # \dontrun{
 library(tidyverse)
 # example tree list
 tl <- dplyr::tibble(
     treeID = c(1:21)
     , tree_x = rnorm(n=21, mean = 458064, sd = 11)
     , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
   )
 # call the function
 tl_lf <- trees_landfire_cbd(tree_list = tl, crs = "32613")
 # what?
 tl_lf %>% class()
 # a list, but what is in it?
 tl_lf %>% names()
 # what's in the trees data?
 tl_lf$tree_list %>% dplyr::glimpse()
 # plot the tree_list spatial points
 tl_lf$tree_list %>% ggplot2::ggplot() +
  ggplot2::geom_sf(ggplot2::aes(color=landfire_cell_kg_per_m3))
 # plot the landfire cbd raster
 tl_lf$landfire_rast %>% terra::plot()
 # we can overlay these
 tl_lf$landfire_rast %>%
   terra::as.data.frame(xy = T) %>%
   ggplot2::ggplot() +
   ggplot2::geom_tile(ggplot2::aes(x=x,y=y,fill=kg_per_m3), color = "gray") +
   ggplot2::geom_sf(
     data = tl_lf$tree_list %>% sf::st_transform(terra::crs(tl_lf$landfire_rast))
     , mapping = ggplot2::aes(color=landfire_cell_kg_per_m3)
   ) +
   ggplot2::scale_color_distiller(palette = "Oranges")
 } # }
```
