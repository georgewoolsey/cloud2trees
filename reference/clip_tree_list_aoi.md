# internal functions to work with tree list data within a specified area of interest (AOI) or "domain" in the terminology of really smart fire modellers.

internal functions to work with tree list data within an AOI. For
example, to use the tree list (e.g. as exported by
[`raster2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/raster2trees.md))
within the QUIC-Fire modelling tool some additional columns are required
and one needs to define a study area and align the DTM and tree list
with this domain.

## Usage

``` r
clip_tree_list_aoi(
  tree_list,
  crs,
  study_boundary,
  bbox_aoi = F,
  buffer = 0,
  reproject_epsg = NULL
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
  which may extend beyond the space with trees. This must be an sf class
  object with a single record. If you need to get the trees within
  multiple different AOI's, then
  [`purrr::map()`](https://purrr.tidyverse.org/reference/map.html) this
  function.

- bbox_aoi:

  logical. Should the study_boundary be transformed to a bounding box
  instead of it's original shape for determining the trees within the
  boundary? If set to true, the bounding box is created prior to
  applying the buffer.

- buffer:

  numeric. Buffer to be applied to the study area prior to determining
  trees within the boundary. Units are determined by the horizontal CRS
  settings of the tree_list data or the CRS of the reproject_epsg.

- reproject_epsg:

  numeric. The EPSG code to reproject the data in prior to buffering and
  clipping. Will determine the projection of the output data.
