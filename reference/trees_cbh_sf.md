# Estimate CBH using tree crown polygons and normalized point cloud data

`trees_cbh_sf()` uses the input tree crown polygons (e.g. as exported by
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
trees_cbh_sf(
  trees_poly,
  norm_las = NULL,
  tree_sample_n = NA,
  tree_sample_prop = NA,
  which_cbh = "lowest",
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
  trees_sample = NA,
  ofile = NA
)
```

## Arguments

- trees_poly:

  sf. A `sf` class object with POLYGON geometry (see
  [`sf::st_geometry_type()`](https://r-spatial.github.io/sf/reference/st_geometry_type.html)),
  the program will use the data "as-is" and only require the `treeID`
  and `tree_height_m` columns. Or the path to a single spatial polygon
  file.

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

- trees_sample:

  data.frame. provide your own tree sample list such as one generated
  from `sample_trees_flist()` that includes the `treeID` column. If
  provided, the tree_sample_n,tree_sample_prop will be ignored

- ofile:

  character or logical. if a character value is provided the output will
  be written to the disk as a csv at the location provided. If set to
  TRUE and a file path was used as the input for `trees_poly`, then a
  csv file will be written to the same location with the same name
  prefixed with "cbh\_". Leave as NA to return a data.frame of the trees
  from tree list from `trees_poly` with CBH values added
