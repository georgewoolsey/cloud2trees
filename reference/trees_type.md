# Estimate forest type for a tree list based on location

`trees_type()` uses the input tree list (e.g. as exported by
[`raster2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/raster2trees.md))
with the columns `treeID`, `tree_x`, `tree_y` to attach species
information using USDA Forest Inventory and Analysis (FIA) codes. If a
spatial data frame of points is the input tree list, then the columns
`tree_x`, `tree_y` are not required.

FIA Forest Type Group Code is attached to each tree in the tree list
based on the spatial overlap with the Forest Type Groups of the
Continental United States dataset [Wilson
2023](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d).

The simplified process for attaching forest type group to a tree is:

- Forest type group 30-m raster (Wilson 2023) was aggregated to 90-m to
  make the data more accessible over the entire continental US

- Nearest neighbor imputation is used to fill forest type data if a tree
  falls inside a no-forest cell in the original data

- The FIA forest type group is applied to a tree based on spatial
  overlap

## Usage

``` r
trees_type(
  tree_list,
  crs = NA,
  study_boundary = NA,
  input_foresttype_dir = NULL,
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

  sf. The boundary of the study area to define the area of the regional
  model. If no boundary given, regional model will be built from
  location of trees in the tree list.

- input_foresttype_dir:

  directory where Forest Type Groups data exists. Use
  [`get_foresttype()`](https://georgewoolsey.github.io/cloud2trees/reference/get_foresttype.md)
  first.

- max_search_dist_m:

  number. Maximum search distance (m) to obtain forest type group data
  for trees in `tree_list` that overlap with non-forest data in the
  original Wilson (2023) data. Larger search distances will increase
  processing time and possibly result in memory issues.

## Value

Returns a list of objects: tree_list = spatial data frame of individual
trees; foresttype_rast = raster of forest types in the area.

## References

- [Forest Type Groups of the Continental United
  States](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d)
  Wilson, B.T. (2023). Forest Type Groups of the Continental United
  States.

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
 tl_type <- trees_type(tree_list = tl, crs = "32613")
 # what?
 tl_type %>% class()
 # a list, but what is in it?
 tl_type %>% names()
 # plot the tree_list spatial points
 tl_type$tree_list %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=forest_type_group))
 # plot the foresttype_rast raster
 tl_type$foresttype_rast %>% terra::plot()
 } # }
```
