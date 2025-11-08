# Estimate tree crown biomass for a tree list using LANDFIRE data

`trees_biomass_landfire()` uses the input tree list (e.g. as exported by
[`raster2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/raster2trees.md))
with the columns `treeID`, `tree_x`, `tree_y` to attach an estimate of
tree crown biomass using LANDFIRE's Forest Canopy Bulk Density (CBD)
data produced jointly by the U.S. Department of Agriculture and U.S.
Department of the Interior.

If a spatial data frame of points is the input tree list, then the
columns `tree_x`, `tree_y` are not required. Other required columns
include:

- `crown_area_m2`, `tree_height_m` (e.g. as exported by
  [`raster2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/raster2trees.md))

- `tree_cbh_m` (e.g. as exported by
  [`trees_cbh()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_cbh.md))

- and one of `dbh_cm`, `dbh_m`, or `basal_area_m2` (e.g. as exported by
  [`trees_dbh()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_dbh.md))

LANDFIRE's Forest Canopy Bulk Density (CBD) data is attached to each
tree in the tree list based on the spatial overlap with the raster data
set (see references). Canopy Bulk Density is mass of flammable material
per unit volume of the tree crown typically expressed in units of mass
per unit volume (e.g., kilograms per cubic meter).

The process for estimating tree crown biomass in kilograms is:

- Nearest neighbor imputation is used to fill LANDFIRE data if a tree
  falls inside a non-forest cell in the original data

- The LANDFIRE estimate of CBD is distributed across the individual
  trees that fall in a raster cell by:

  1.  At the stand level (i.e. raster cell), aggregate the tree level
      data within the stand to obtain:

  - `mean_crown_length_m = mean(crown_length_m)`, where tree
    `crown_length_m = tree_height_m - tree_cbh_m`

  - `sum_crown_volume_m3 = sum(crown_volume_m3)`, where tree
    `crown_volume_m3 = (4/3) * pi * ((crown_length_m/2)) * ((crown_dia_m/2)^2)`

  1.  At the stand level (i.e. raster cell), determine the area of the
      stand that overlaps (`overlap_area_m2`) with the AOI defined as
      the `study_boundary` parameter (see below) or the bounding box of
      all the trees

  2.  At the stand level (i.e. raster cell), get the LANDFIRE estimate
      of CBD in kilograms per cubic meter (`landfire_stand_kg_per_m3`)

  3.  At the stand level (i.e. raster cell), get canopy fuel loading
      (CFL) in kilograms per square meter
      (`kg_per_m2 = mean_crown_length_m * landfire_stand_kg_per_m3`)

  4.  At the stand level (i.e. raster cell), get the stand biomass in
      kilograms (`biomass_kg = kg_per_m2 * overlap_area_m2`)

  5.  At the stand level (i.e. raster cell), the single tree CBD in
      kilograms per cubic meter will be a constant
      (`landfire_tree_kg_per_m3 = biomass_kg / sum_crown_volume_m3`)

  6.  Attach the the single tree CBD in kilograms per cubic meter to the
      tree level based on raster cell spatial overlap

  7.  Calculate individual tree crown mass in kilograms as
      `landfire_crown_biomass_kg = landfire_tree_kg_per_m3 * crown_volume_m3`

## Usage

``` r
trees_biomass_landfire(
  tree_list,
  crs = NA,
  study_boundary = NA,
  input_landfire_dir = NULL,
  max_crown_kg_per_m3 = 2
)
```

## Arguments

- tree_list:

  data.frame. A data frame with the columns `treeID`, `tree_x`, and
  `tree_y`. If an `sf` class object with POINT geometry (see
  [`sf::st_geometry_type()`](https://r-spatial.github.io/sf/reference/st_geometry_type.html)),
  the program will use the data "as-is" and only require the `treeID`
  column. Other required columns include:

  - `crown_area_m2`, `tree_height_m` (e.g. as exported by
    [`raster2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/raster2trees.md))

  - `tree_cbh_m` (e.g. as exported by
    [`trees_cbh()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_cbh.md))

  - and one of `dbh_cm`, `dbh_m`, or `basal_area_m2` (e.g. as exported
    by
    [`trees_dbh()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_dbh.md))

- crs:

  string. A crs string as returned from
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)
  or the EPSG code of the x,y coordinates. Defaults to the crs of the
  `tree_list` data if of class "sf".

- study_boundary:

  sf. The boundary of the study area to define the area of interest
  which may extend beyond the space with trees. If no boundary given,
  the AOI will be built from location of trees in the tree list.

- input_landfire_dir:

  directory where LANDFIRE CBD data exists. Use
  [`get_landfire()`](https://georgewoolsey.github.io/cloud2trees/reference/get_landfire.md)
  first.

- max_crown_kg_per_m3:

  numeric. the maximum CBD of the tree crown in kilograms per cubic
  meter. Values above this limit will be set at the median value for the
  area using only stands that have CBD values lower than this limit. The
  default value of 2 kilograms per cubic meter was based on [Mell et al.
  (2009)](https://doi.org/10.1016/j.combustflame.2009.06.015) who found
  the dry bulk density of the tree crown was 2.6 kilograms per cubed
  meter using Douglas-fir trees grown on Christmas tree farms. Set this
  parameter to a large value (e.g. 1e10) or NULL to avoid limiting tree
  crown CBD.

## Value

Returns a list of objects: tree_list = spatial data frame of individual
trees ; stand_cell_data = data frame of stands/cells in same projection
as the LANDFIRE raster data

See code in examples.

## References

- [LANDFIRE Forest Canopy Bulk Density
  (CBD)](https://landfire.gov/fuel/cbd) U.S. Department of Agriculture
  and U.S. Department of the Interior.

- <https://doi.org/10.1016/j.combustflame.2009.06.015> Mell, W.,
  Maranghides, A., McDermott, R., & Manzello, S. L. (2009). Numerical
  simulation and experiments of burning douglas fir trees. Combustion
  and Flame, 156(10), 2023-2041.

## Examples

``` r
 if (FALSE) { # \dontrun{
library(tidyverse)
library(sf)
my_n <- 122
# fake tree list
tl <- dplyr::tibble(
    treeID = c(1:my_n)
    , tree_x = rnorm(n=my_n, mean = 458064, sd = 33)
    , tree_y = rnorm(n=my_n, mean = 4450074, sd = 33)
    , tree_height_m = rnorm(n=my_n, mean = 10, sd = 7)
  ) %>%
  dplyr::mutate(
    tree_height_m = ifelse(tree_height_m<1.37, 1.37, tree_height_m) # above DBH
    , crown_area_m2 = 0.47+(0.49*tree_height_m)
    , tree_cbh_m = 0.72+(0.53*tree_height_m)
    , dbh_cm = -2.3+(2.14*tree_height_m)
  )
# how does our fake tree list look?
tl %>% dplyr::glimpse()
tl %>% ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(x = dbh_cm, y = tree_height_m))
# call the function
tl_landfire <- trees_biomass_landfire(tree_list = tl, crs = "32613")
# what is in it?
tl_landfire %>% names()
# look at the trees
tl_landfire$tree_list %>% dplyr::glimpse()
# look at the stand
tl_landfire$stand_cell_data %>% dplyr::filter(!is.na(trees)) %>% dplyr::glimpse()
# the stand CBD
 tl_landfire$stand_cell_data %>%
   dplyr::filter(trees>0) %>%
   sf::st_drop_geometry() %>%
   dplyr::count(landfire_stand_kg_per_m3)
# get the projection for the stand cell data
epsg_code <- tl_landfire$stand_cell_data$rast_epsg_code[1] %>% as.numeric()
# plot the stand cell data with trees overlaid
tl_landfire$stand_cell_data %>%
  ggplot2::ggplot() +
  ggplot2::geom_tile(ggplot2::aes(x=x,y=y,fill = landfire_stand_kg_per_m3), color = "gray44") +
  ggplot2::geom_text(ggplot2::aes(x=x,y=y,label = trees), color = "white") +
  ggplot2::geom_sf(
    data = tl_landfire$tree_list %>% sf::st_transform(crs = epsg_code)
    , ggplot2::aes(color = landfire_crown_biomass_kg)
  ) +
  ggplot2::labs(fill="stand kg/m3", color = "tree crown kg", caption = "# trees shown in cell") +
  ggplot2::scale_fill_viridis_c(option = "rocket", na.value = "gray", direction = -1) +
  ggplot2::scale_color_viridis_c(option = "viridis", na.value = "gray22", begin = 0.6) +
  ggplot2::theme_void()
 } # }
```
