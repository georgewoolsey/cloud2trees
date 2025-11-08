# Use outputs from [`cloud2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/cloud2trees.md) to generate inputs for LANL TREES program

`cloud2trees_to_lanl_trees()` uses the output from
[`cloud2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/cloud2trees.md)
to generate inputs for LANL TREES program as a pathway to fire modeling
with Quic-Fire

The primary input is a directory with outputs from
[`cloud2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/cloud2trees.md).
The default directory written by
[`cloud2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/cloud2trees.md)
is `point_cloud_processing_delivery` which must contain (at a minimum):

- DTM raster with a name formatted as: "dtm_xx.tif"

- Tree list data with the name formatted as :
  "final_detected_tree_tops.gpkg" (tree points) or
  "final_detected_crowns.gpkg" (tree crowns)

- A study area spatial file that can be read with the `sf` package (see
  [`sf::st_drivers()`](https://r-spatial.github.io/sf/reference/st_drivers.html))

## Usage

``` r
cloud2trees_to_lanl_trees(
  input_dir,
  study_boundary = NA,
  bbox_aoi = T,
  buffer = 0,
  topofile = "flat",
  cbd_method = "landfire",
  output_dir = tempdir(),
  fuel_litter = list(ilitter = 0, lrho = 4.667, lmoisture = 0.06, lss = 5e-04, ldepth =
    0.06),
  fuel_grass = list(igrass = 0, grho = 1.17, gmoisture = 0.06, gss = 5e-04, gdepth =
    0.27)
)
```

## Arguments

- input_dir:

  directory with outputs from
  [`cloud2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/cloud2trees.md).
  The default directory written by
  [`cloud2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/cloud2trees.md)
  is `point_cloud_processing_delivery`

- study_boundary:

  sf. The boundary of the study area which is used to determine the
  outputs

- bbox_aoi:

  logical. Should the study_boundary be transformed to a bounding box
  instead of it's original shape for determining the objects within the
  boundary?

- buffer:

  numeric. Buffer to be applied to the study area prior to determining
  objects within the boundary. Units are determined by the horizontal
  CRS settings of the tree list data

- topofile:

  character. one of `"flat"` or `"dtm"`:

  - "flat" - always flat for QUIC-Fire

  - "dtm" - uses the topo.dat file created based on the DTM; potentially
    for FIRETEC

- cbd_method:

  character. one of `"landfire"` or `"cruz"`:

  - Tree crown biomass method:

    - "landfire" - based on LANDFIRE's Forest Canopy Bulk Density (CBD)
      data
      ([`trees_biomass_landfire()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_biomass_landfire.md))

    - "cruz" - based on Cruz et al. (2003) canopy fuel stratum equations
      ([`trees_biomass_cruz()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_biomass_cruz.md))

- output_dir:

  parent directory where new folder `lanl_trees_delivery` will be
  written for exports

- fuel_litter:

  list. a [`list()`](https://rdrr.io/r/base/list.html) or numeric vector
  [`c()`](https://rdrr.io/r/base/c.html). see default.

  - must have parameters in order:

    - "ilitter" : 0 = no litter, 1 = litter

    - "lrho" : litter bulk density (kg/m3)

    - "lmoisture : litter moisture (percent on 0-1 scale)

    - "lss" : litter sizescale (m)

    - "ldepth" : litter depth (m)

- fuel_grass:

  list. a [`list()`](https://rdrr.io/r/base/list.html) or numeric vector
  [`c()`](https://rdrr.io/r/base/c.html). see default.

  - must have parameters in order:

    - "igrass" : 0 = no grass, 1 = grass

    - "grho" : grass bulk density (kg/m3)

    - "gmoisture" : grass moisture (percent on 0-1 scale)

    - "gss" : grass sizescale (m)

    - "gdepth" : grass depth (m)

## Value

Returns a list of objects:

- "tree_list" = the cropped tree list based on the study area extent
  with customized settings

- "aoi" = the study area extent with customized settings

- "dtm" = the cropped DTM based on the study area extent with customized
  settings

- "domain_path" = the path to the "Lidar_Bounds.geojson" file

- "topofile_path" = the path to the "topo.dat" file

- "fuellist_path" = the path to the TREES program "fuellist" file

- "treelist_path" = the path to the "Cloud2Trees_TreeList.txt" file

## References

- <https://github.com/lanl/Trees/>

- [Quic-Fire](https://doi.org/10.1016/j.envsoft.2019.104616)

## Examples

``` r
 if (FALSE) { # \dontrun{
 # test las file but this could also be a directory path with >1 .las|.laz files
 i <- system.file("extdata", "MixedConifer.laz", package="lidR")
 # set the dir to write the output
 my_dir <- tempdir()
 # run it
 cloud2trees_ans <- cloud2trees::cloud2trees(
   output_dir = my_dir
   , input_las_dir = i
   # turn on all of the attribute estimations
   , estimate_tree_dbh = T
   , estimate_tree_type = T
   , estimate_tree_hmd = T
   , hmd_tree_sample_prop = 0.5
   , hmd_estimate_missing_hmd = T
   , estimate_biomass_method = "landfire"
   , estimate_tree_cbh = T
   , cbh_tree_sample_prop = 0.3
   , cbh_estimate_missing_cbh = T
 )
 # generate a fake study_boundary we know overlaps the tree list
 my_aoi <- cloud2trees_ans$treetops_sf %>%
   sf::st_bbox() %>%
   sf::st_as_sfc() %>%
   sf::st_buffer(-10) %>%
   sf::st_as_sf()
 # plot it
 ggplot2::ggplot() +
   ggplot2::geom_sf(data = cloud2trees_ans$treetops_sf) +
   ggplot2::geom_sf(data = my_aoi, fill = NA, color = "blue")
 # cloud2trees::cloud2trees() wrote the `point_cloud_processing_delivery` folder
 cloud2trees_output_dir <- file.path(my_dir,"point_cloud_processing_delivery")
 list.files(cloud2trees_output_dir)
 list.dirs(cloud2trees_output_dir, recursive = F)
 # now cloud2trees_to_lanl_trees()
 cloud2trees_to_lanl_trees_ans <- cloud2trees_to_lanl_trees(
   input_dir = cloud2trees_output_dir
   , study_boundary = my_aoi
   , bbox_aoi = F
   , buffer = 0
   , topofile = "flat"
   , cbd_method = "landfire"
   , output_dir = cloud2trees_output_dir
 )
 cloud2trees_to_lanl_trees_ans %>% names()
 list.dirs(cloud2trees_output_dir, recursive = F)
 lanl_trees_output_dir <- file.path(cloud2trees_output_dir,"lanl_trees_delivery")
 list.files(lanl_trees_output_dir)
 } # }
```
