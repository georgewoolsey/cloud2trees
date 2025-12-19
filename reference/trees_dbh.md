# Estimate DBH for a tree list based on height

`trees_dbh()` uses the input tree list (e.g. as exported by
[`raster2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/raster2trees.md))
with the columns `treeID`, `tree_x`, `tree_y`, and `tree_height_m` to
estimate tree DBH.

A regional model of height estimating DBH is determined by the process:

- Use the TreeMap
  ([`get_treemap()`](https://georgewoolsey.github.io/cloud2trees/reference/get_treemap.md))
  FIA plot data in the area of the tree list to estimate the height-DBH
  allometry relationship

- Use the height predicting DBH model built from the FIA data to predict
  DBH based on tree height in the tree list

If training data is provided in `treels_dbh_locations` as returned from
[`treels_stem_dbh()`](https://georgewoolsey.github.io/cloud2trees/reference/treels_stem_dbh.md):

- The regional model using FIA plot data is used to filter the DBH
  training data estimated from the point cloud

- The training data is used to estimate the height-DBH allometry
  relationship

- Use the height predicting DBH model built from the point cloud
  training data to predict DBH based on tree height in the tree list

## Usage

``` r
trees_dbh(
  tree_list,
  crs = NA,
  study_boundary = NA,
  dbh_model_regional = "cr",
  dbh_model_local = "lin",
  treels_dbh_locations = NA,
  boundary_buffer = 50,
  input_treemap_dir = NULL,
  outfolder = tempdir()
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

- study_boundary:

  sf. The boundary of the study are to define the area of the regional
  model. If no boundary given, regional model will be built from
  location of trees in the tree list.

- dbh_model_regional:

  string. Set the model to use for regional dbh-height allometry based
  on FIA tree measurements. Can be "cr" for the Chapman-Richards formula
  (default) or "power" for power function

- dbh_model_local:

  string. Set the model to use for local dbh-height allometry based on
  provided DBH training data in `treels_dbh_locations`. Can be "rf" for
  random forest or "lin" for linear

- treels_dbh_locations:

  sf. Return from
  [`treels_stem_dbh()`](https://georgewoolsey.github.io/cloud2trees/reference/treels_stem_dbh.md).
  Must also provide crown polygons (as returned from
  [`raster2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/raster2trees.md))
  in the `tree_list` data as an `sf` class object with POLYGON geometry
  (see
  [`sf::st_geometry_type()`](https://r-spatial.github.io/sf/reference/st_geometry_type.html))
  If a valid file is provided, will make DBH predictions based on this
  training data instead of from the regional model from the FIA data

- boundary_buffer:

  numeric. Set the buffer (m) for the study area boundary to filter the
  FIA plot data based on TreeMap

- input_treemap_dir:

  directory where Treemap 2016 exists. Use
  [`get_treemap()`](https://georgewoolsey.github.io/cloud2trees/reference/get_treemap.md)
  first.

- outfolder:

  string. The path of a folder to write the model data to

- dbh_model:

  **\[deprecated\]** Use the `dbh_model_regional` or `dbh_model_local`
  argument instead.

## Value

Returns a spatial data frame of individual trees.

## References

- <https://doi.org/10.2737/RDS-2025-0032> Houtman, Rachel M.;
  Leatherman, Lila S. T.; Zimmer, Scott N.; Housman, Ian W.; Shrestha,
  Abhinav; Shaw, John D.; Riley, Karin L. 2025. TreeMap 2022 CONUS: A
  tree-level model of the forests of the conterminous United States
  circa 2022. Fort Collins, CO: Forest Service Research Data Archive.

- <https://doi.org/10.3390/f13122077> Tinkham et al. (2022). Modeling
  the missing DBHs: Influence of model form on UAV DBH characterization.
  Forests, 13(12), 2077.

## Examples

``` r
 if (FALSE) { # \dontrun{
 library(tidyverse)
 # example tree list
 tl <- dplyr::tibble(
     treeID = c(1:21)
     , tree_x = rnorm(n=21, mean = 458064, sd = 11)
     , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
     , tree_height_m = exp(rgamma(n = 21, shape = (7/4)^2, rate = (4^2)/7))
   )
 # save our output somewhere (not required)
 outdir <- tempdir()
 # call the function
 tl_dbh <- trees_dbh(tree_list = tl, crs = "32613", outfolder = outdir)
 # what?
 tl_dbh %>% class()
 tl_dbh %>% dplyr::select(tidyselect::contains("dbh_cm")) %>% dplyr::glimpse()
 tl_dbh %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=dbh_cm))
 # what outputs did we get?
 list.files(outdir)
 # cloud2trees::trees_dbh() saved the FIA-measured trees used to train the allometric model
 read.csv(file.path(outdir, "regional_dbh_height_model_training_data.csv")) %>%
   summary()
 # cloud2trees::trees_dbh() saved the actual allometric model
 # let's load and review
 dbh_mod_temp <- readRDS(file.path(outdir, "regional_dbh_height_model.rds"))
 # what is this?
 dbh_mod_temp %>% class()
 # we can draw fit curves with probability bands using the tidybayes package
 library(tidybayes)
 # define our height range to predict over
 dplyr::tibble(tree_height_m = seq(from = 0, to = 25, by = 1)) %>%
   tidybayes::add_epred_draws(dbh_mod_temp, ndraws = 2000) %>%
   ggplot2::ggplot(ggplot2::aes(x = tree_height_m)) +
     tidybayes::stat_lineribbon(
       ggplot2::aes(y = .epred, color = "estimate")
       , .width = c(0.5,0.95)
       , lwd = 0.6
     ) +
     ggplot2::scale_fill_brewer(palette = "Oranges") +
     ggplot2::scale_color_manual(values = c("gray33")) +
     ggplot2::labs(x = "tree ht. (m)", y = "est. tree DBH (cm)", color = "") +
     ggplot2::scale_x_continuous(limits = c(0,NA), breaks = scales::extended_breaks(n=11)) +
     ggplot2::scale_y_continuous(limits = c(0,NA), breaks = scales::extended_breaks(n=11)) +
     ggplot2::theme_light()
 } # }
```
