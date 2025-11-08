# Estimate tree biomass (or crown biomass) for a tree list

`trees_biomass()` streamlines the process for estimating individual tree
biomass in kilograms, or the component biomass of the tree crown in
kilograms. Users can select one or all of the following methods
available in the package for estimating biomass:

- Tree crown biomass in kilograms:

  - "landfire" - based on LANDFIRE's Forest Canopy Bulk Density (CBD)
    data
    ([`trees_biomass_landfire()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_biomass_landfire.md))

  - "cruz" - based on Cruz et al. (2003) canopy fuel stratum equations
    ([`trees_biomass_cruz()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_biomass_cruz.md))

- Tree total above ground biomass in kilograms:

  - coming soon

If multiple methods are selected (e.g. `method = c("cruz","landfire")`),
then the program will compile the biomass estimates and return one tree
list.

## Usage

``` r
trees_biomass(
  tree_list,
  crs = NA,
  study_boundary = NA,
  input_landfire_dir = NULL,
  input_foresttype_dir = NULL,
  method = NA,
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

- input_foresttype_dir:

  directory where Forest Type Groups data exists. Use
  [`get_foresttype()`](https://georgewoolsey.github.io/cloud2trees/reference/get_foresttype.md)
  first.

- method:

  character. one (e.g. `"landfire"`) or multiple (e.g.
  `c("cruz","landfire")`) of the following biomass estimation methods:

  - Tree crown biomass in kilograms:

    - "landfire" - based on LANDFIRE's Forest Canopy Bulk Density (CBD)
      data
      ([`trees_biomass_landfire()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_biomass_landfire.md))

    - "cruz" - based on Cruz et al. (2003) canopy fuel stratum equations
      ([`trees_biomass_cruz()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_biomass_cruz.md))

  - Tree total above ground biomass in kilograms:

    - coming soon

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
trees ; stand_cell_data_landfire = data frame of stands/cells in same
projection as the LANDFIRE raster data ; stand_cell_data_cruz = data
frame of stands/cells in same projection as the FIA forest type group
raster data

See code in examples.

## References

See references in:

- [`trees_biomass_landfire()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_biomass_landfire.md)

- [`trees_biomass_cruz()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_biomass_cruz.md)

## Examples

``` r
 if (FALSE) { # \dontrun{
library(tidyverse)
library(sf)
# use the tree list that ships with the package
f <- system.file(package = "cloud2trees", "extdata", "crowns_poly.gpkg")
tl <- sf::st_read(f)
tl %>% dplyr::glimpse()
# call trees_biomass and get multiple biomass estimates
trees_biomass_ans <- trees_biomass(tree_list = tl, method = c("landfire","cruz"))
# what did we get back?
trees_biomass_ans %>% names()
# check out the tree list
trees_biomass_ans$tree_list %>% dplyr::glimpse()
# check out the landfire stand data
trees_biomass_ans$stand_cell_data_landfire %>% dplyr::filter(trees>0) %>% dplyr::glimpse()
# plot tree landfire crown biomass estimate
trees_biomass_ans$tree_list %>%
  ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = tree_height_m
      , y = landfire_crown_biomass_kg
      , color = crown_area_m2
    )
  ) +
  ggplot2::geom_point()
# plot tree cruz crown biomass estimate
trees_biomass_ans$tree_list %>%
  ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = tree_height_m
      , y = cruz_crown_biomass_kg
      , color = crown_area_m2
    )
  ) +
  ggplot2::geom_point()
# plot tree landfire vs. cruz crown biomass estimate
trees_biomass_ans$tree_list %>%
  ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = landfire_crown_biomass_kg, y = cruz_crown_biomass_kg
    )
  ) +
  ggplot2::geom_abline(lwd = 1.5) +
  ggplot2::geom_smooth(method = "lm", se=F, color = "gray", linetype = "dashed") +
  ggplot2::geom_point(ggplot2::aes(color = tree_height_m)) +
  ggplot2::scale_x_continuous(
    limits = c(0
      , max(trees_biomass_ans$tree_list$cruz_crown_biomass_kg)
    )
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0
      , max(trees_biomass_ans$tree_list$cruz_crown_biomass_kg)
    )
  )
# get the projection for the stand cell data
epsg_code <- trees_biomass_ans$stand_cell_data_landfire$rast_epsg_code[1] %>% as.numeric()
# plot the stand cell data with trees overlaid
trees_biomass_ans$stand_cell_data_landfire %>%
  dplyr::filter(trees>0) %>%
  ggplot2::ggplot() +
  ggplot2::geom_tile(ggplot2::aes(x=x,y=y,fill = landfire_stand_kg_per_m3), color = "gray44") +
  ggplot2::geom_text(ggplot2::aes(x=x,y=y,label = trees), color = "white") +
  ggplot2::geom_sf(
    data = trees_biomass_ans$tree_list %>% sf::st_transform(crs = epsg_code)
    , ggplot2::aes(color = cruz_crown_biomass_kg)
  ) +
  ggplot2::labs(
    fill="stand kg/m3", color = "landfire\ncrown kg"
    , caption = "# trees shown in cell"
  ) +
  ggplot2::scale_fill_viridis_c(option = "rocket", na.value = "gray", direction = -1) +
  ggplot2::scale_color_viridis_c(option = "viridis", na.value = "gray22", begin = 0.6) +
  ggplot2::theme_void()
 } # }
```
