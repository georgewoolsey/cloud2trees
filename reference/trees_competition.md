# Calculate competition metrics for a tree list

`trees_competition()` uses the input tree list (e.g. as exported by
[`raster2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/raster2trees.md))
with the columns `treeID`, `tree_x`, `tree_y`, and `tree_height_m` to
calculate competition metrics at the tree level.

Competition metrics returned include:

- Distance to the nearest neighbor (`comp_dist_to_nearest_m`)

- Trees per ha within a 5m radius (`comp_trees_per_ha`)

- The relative tree height (`comp_relative_tree_height` = height_tree /
  height_max) within a 5m radius where a value of 1 indicates the
  tallest tree.

## Usage

``` r
trees_competition(
  tree_list,
  crs = NA,
  competition_buffer_m = 5,
  study_boundary = NA,
  search_dist_max = 10
)
```

## Arguments

- tree_list:

  data.frame. A data frame with the columns `treeID`, `tree_x`,
  `tree_y`, and `tree_height_m`. If an `sf` class object with POINT
  geometry (see
  [`sf::st_geometry_type()`](https://r-spatial.github.io/sf/reference/st_geometry_type.html)),
  the program will use the data "as-is" and only require the `treeID`
  and `tree_height_m` columns.

- crs:

  string. A crs string as returned from
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)
  or the EPSG code of the x,y coordinates. Defaults to the crs of the
  `tree_list` data if of class "sf".

- competition_buffer_m:

  number. Set buffer around tree (m) to calculate competition metrics

- study_boundary:

  sf. If you want to scale per ha calculations, provide the geography of
  the study boundary

- search_dist_max:

  number. Maximum search distance (m) to nearest tree. Larger search
  distances will increase processing time and possibly result in memory
  issues. If no competition trees are found within this distance, the
  return column `comp_dist_to_nearest_m` = `search_dist_max` parameter.

## Value

Returns a spatial data frame of individual trees.

## References

<https://doi.org/10.3390/f13122077> Tinkham et al. (2022). Modeling the
missing DBHs: Influence of model form on UAV DBH characterization.
Forests, 13(12), 2077.

## Examples

``` r
 if (FALSE) { # \dontrun{
 # example tree list
 tl <- dplyr::tibble(
     treeID = c(1:21)
     , tree_x = rnorm(n=21, mean = 458064, sd = 11)
     , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
     , tree_height_m = exp(rgamma(n = 21, shape = (7/4)^2, rate = (4^2)/7))
   )
 # call the function
 tl_comp <- trees_competition(tree_list = tl, crs = "32613")
 # what?
 tl_comp %>% class()
 tl_comp %>% dplyr::select(tidyselect::starts_with("comp_")) %>% dplyr::glimpse()
 tl_comp %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=comp_dist_to_nearest_m))
 } # }
```
