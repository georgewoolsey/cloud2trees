# Estimate CBH using tree crown polygons and normalized point cloud data

`trees_cbh()` uses the input tree crown polygons (e.g. as exported by
[`raster2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/raster2trees.md))
with the columns `treeID` and `tree_height_m` to estimate tree CBH using
height normalized point cloud data (e.g. as exported by
[`cloud2raster()`](https://georgewoolsey.github.io/cloud2trees/reference/cloud2raster.md)).

CBH is extracted directly from the height normalized point cloud using
the process outlined in Viedma et al. (2024) and implemented via
[`ladderfuelsr_cbh()`](https://georgewoolsey.github.io/cloud2trees/reference/ladderfuelsr_cbh.md).

There are likely to be trees for which there is insufficient data in the
point cloud to successfully estimate CBH. The user can elect to estimate
missing CBH values which is accomplished via:

- Attempt to extract CBH from the sample of trees elected by the user
  (`tree_sample_n`,`tree_sample_prop` parameter) using
  [`ladderfuelsr_cbh()`](https://georgewoolsey.github.io/cloud2trees/reference/ladderfuelsr_cbh.md)

- Successfully extracted CBH trees become training data used to estimate
  the height-CBH allometry relationship that is spatially informed using
  the relative tree location compared to the training data

- The height and location predicting CBH model built from the point
  cloud training data is used to predict CBH for the non-training (i.e.
  missing CBH) data

## Usage

``` r
trees_cbh(
  trees_poly,
  norm_las = NULL,
  tree_sample_n = NA,
  tree_sample_prop = NA,
  which_cbh = "lowest",
  estimate_missing_cbh = TRUE,
  min_vhp_n = 3,
  voxel_grain_size_m = 1,
  dist_btwn_bins_m = 1,
  min_fuel_layer_ht_m = 1,
  lad_pct_gap = 25,
  lad_pct_base = 25,
  num_jump_steps = 1,
  min_lad_pct = 10,
  frst_layer_min_ht_m = 1,
  force_same_crs = F,
  outfolder = tempdir()
)
```

## Arguments

- trees_poly:

  must be one of the following that has required attributes `treeID` and
  `tree_height_m`:

  - `sf` class object with POLYGON geometry (see
    [`sf::st_geometry_type()`](https://r-spatial.github.io/sf/reference/st_geometry_type.html)).
    Recommended for smaller tree lists (e.g. \<100k) that can fit in
    memory.

  - character vector with the path to a single or multiple spatial files
    that can be read by
    [`sf::st_read()`](https://r-spatial.github.io/sf/reference/st_read.html)
    and have POLYGON geometry. Recommended for large tree lists (e.g.
    100k+) that might cause memory issues.

  - character with the path to a directory that has
    "final_detected_crowns\*" files from
    [`cloud2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/cloud2trees.md)
    or
    [`raster2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/raster2trees.md).
    Recommended for large tree lists (e.g. 100k+) that might cause
    memory issues.

- norm_las:

  character. a directory with nomalized las files, the path of a single
  .laz\|.las file", -or- an object of class `LAS`. It is your
  responsibility to ensure that the point cloud is projected the same as
  the `trees_poly` data

- tree_sample_n, tree_sample_prop:

  numeric. Provide either `tree_sample_n`, the number of trees, or
  `tree_sample_prop`, the proportion of the trees to attempt to extract
  a CBH from the point cloud for. If neither are supplied,
  `tree_sample_n = 333` will be used. If both are supplied,
  `tree_sample_n` will be used. Increasing `tree_sample_prop` toward
  one (1) will increase the processing time, perhaps significantly
  depending on the number of trees in the `trees_poly` data.

- which_cbh:

  character. One of: "lowest"; "highest"; or "max_lad". See Viedma et
  al. (2024) reference.

  - "lowest" - Height of the CBH of the segmented tree based on the last
    distance found in its profile

  - "highest" - Height of the CBH of the segmented tree based on the
    maximum distance found in its profile

  - "max_lad" - Height of the CBH of the segmented tree based on the
    maximum LAD percentage

- estimate_missing_cbh:

  logical. even if the `tree_sample_prop` parameter is set to "1", it is
  not likely that CBH will be extracted successfully from every tree.
  Should the missing CBH values be estimated using the tree height and
  location information based on trees for which CBH is successfully
  extracted?

- min_vhp_n:

  numeric. the minimum number of vertical height profiles (VHPs) needed
  to estimate a CBH.

- voxel_grain_size_m:

  numeric. horizontal resolution (suggested 1 meter for lad profiles and
  10 meters for LAI maps). See `grain.size` in
  [`leafR::lad.voxels()`](https://rdrr.io/pkg/leafR/man/lad.voxels.html)

- dist_btwn_bins_m:

  numeric. value for the actual height bin step (in meters). See `step`
  in
  [`LadderFuelsR::get_gaps_fbhs()`](https://rdrr.io/pkg/LadderFuelsR/man/get_gaps_fbhs.html)

- min_fuel_layer_ht_m:

  numeric. value for the actual minimum base height (in meters). See
  `min_height` in
  [`LadderFuelsR::get_gaps_fbhs()`](https://rdrr.io/pkg/LadderFuelsR/man/get_gaps_fbhs.html)

- lad_pct_gap:

  numeric. value of the percentile threshold used to identify gaps
  (default percentile 25th). See `perc_gap` in
  [`LadderFuelsR::get_gaps_fbhs()`](https://rdrr.io/pkg/LadderFuelsR/man/get_gaps_fbhs.html)

- lad_pct_base:

  numeric. value of the percentile threshold used to identify fuels
  layers base height (default percentile 25th). See `perc_base` in
  [`LadderFuelsR::get_gaps_fbhs()`](https://rdrr.io/pkg/LadderFuelsR/man/get_gaps_fbhs.html)

- num_jump_steps:

  numeric. value for the number of height bin steps that can be jumped
  to reshape fuels layers. See `number_steps` in
  [`LadderFuelsR::get_real_fbh()`](https://rdrr.io/pkg/LadderFuelsR/man/get_real_fbh.html)

- min_lad_pct:

  numeric. value for the minimum required LAD percentage in a fuel
  layer. See `threshold` in
  [`LadderFuelsR::get_layers_lad()`](https://rdrr.io/pkg/LadderFuelsR/man/get_layers_lad.html)

- frst_layer_min_ht_m:

  numeric. value for the depth height of the first fuel layer. If the
  first fuel layer has the maximum LAD and its depth is greater than the
  indicated value, then this fuel layer is considered as the CBH of the
  tree. On the contrary, if its depth is \<= the value, the CBH with
  maximum LAD will be the second fuel layer, although it has not the
  maximum LAD. See `hdepth1_height` in
  [`LadderFuelsR::get_cbh_metrics()`](https://rdrr.io/pkg/LadderFuelsR/man/get_cbh_metrics.html)

- force_same_crs:

  logical. force the same crs between the point cloud and polygon if
  confident that data are in same projection. data created by a
  `cloud2trees` pipeline (e.g.
  [`cloud2raster()`](https://georgewoolsey.github.io/cloud2trees/reference/cloud2raster.md))
  will always have the same projection even if not recognized by `lidR`
  functions

- outfolder:

  string. The path of a folder to write the model data to. Note, in the
  actual missing value estimation many RF models are estimated and model
  averaging is used. However, only the first estimated model is saved in
  this export which does not fully represent the process used to fill in
  missing values.

## Value

Returns a spatial data frame of individual trees.

## References

- <https://doi.org/10.1111/2041-210X.14427> Viedma, O., Silva, C. A.,
  Moreno, J. M., & Hudak, A. T. (2024). LadderFuelsR: A new automated
  tool for vertical fuel continuity analysis and crown base height
  detection using light detection and ranging. Methods in Ecology and
  Evolution. <https://github.com/olgaviedma/LadderFuelsR>

- <https://doi.org/10.3390/rs11010092> Almeida, D. R. A. D., Stark, S.
  C., Shao, G., Schietti, J., Nelson, B. W., Silva, C. A., ... &
  Brancalion, P. H. S. (2019). Optimizing the remote detection of
  tropical rainforest structure with airborne lidar: Leaf area profile
  sensitivity to pulse density and spatial sampling. Remote Sensing,
  11(1), 92. <https://github.com/DRAAlmeida/leafR>

## Examples

``` r
 if (FALSE) { # \dontrun{
 library(tidyverse)
 library(sf)
 # example tree crown polygons
 f <- system.file(package = "cloud2trees","extdata","crowns_poly.gpkg")
 crowns <- sf::st_read(f, quiet = T)
 # example normalized las files are in this directory
 norm_d <- system.file(package = "cloud2trees","extdata","norm_las")
 # now run the trees_cbh()
 trees_cbh_ans <- trees_cbh(
    trees_poly = crowns
    , norm_las = norm_d
    , tree_sample_n = 44
    , estimate_missing_cbh = T
    , force_same_crs = T
   )
 # what?
 trees_cbh_ans %>% class()
 trees_cbh_ans %>% dplyr::select(treeID,tidyselect::contains("cbh")) %>% dplyr::glimpse()
 # spatial polygons
 trees_cbh_ans %>% ggplot2::ggplot() +
   ggplot2::geom_sf(ggplot2::aes(fill=tree_cbh_m,color=is_training_cbh))
 # relationship between height and cbh
 trees_cbh_ans %>%
    ggplot2::ggplot(
      ggplot2::aes(x = tree_height_m, y = tree_cbh_m, color=is_training_cbh)
     ) +
    ggplot2::geom_point()
 # tabulate training data
 trees_cbh_ans %>%
   sf::st_drop_geometry() %>%
   dplyr::count(is_training_cbh)
 #### try a file list
 #### Recommended for large tree lists (e.g. 100k+) that might cause memory issues.
 # we'll split the crowns
 # as is done automatically for tree lists >250k by raster2trees() and cloud2trees()
 crowns <- crowns %>%
   dplyr::mutate(
     # makes 2 groups of data
     grp = ceiling(dplyr::row_number()/(dplyr::n()/2))
   )
 # make file names
 my_dir <- tempdir()
 fnm_1 <- file.path(my_dir, "crowns1.gpkg")
 fnm_2 <- file.path(my_dir, "crowns2.gpkg")
 fnm_1
 # write the data
 sf::st_write(crowns %>% dplyr::filter(grp==1), dsn = fnm_1, append = F) # grp 1
 sf::st_write(crowns %>% dplyr::filter(grp==2), dsn = fnm_2, append = F) # grp 2
 # try trees_cbh with our file list
 flist <- c(fnm_1,fnm_2)
 # now run the trees_cbh()
 trees_cbh_ans2 <- trees_cbh(
  trees_poly = flist
  , norm_las = norm_d
  , tree_sample_n = 44
  , estimate_missing_cbh = T
  , force_same_crs = T
 )
 # tabulate training data
 trees_cbh_ans2 %>%
   sf::st_drop_geometry() %>%
   dplyr::count(is_training_cbh)
 } # }
```
