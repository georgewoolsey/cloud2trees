# Use raw .las\|.laz files to generate CHM, DTM, and a tree list

`cloud2trees()` is an all-in-one function to process raw .las\|.laz
files to generate a CHM raster (.tif), a DTM raster (.tif), and a tree
list with tree location, height, and DBH. The order of operations is:

- Generate a CHM from the point cloud using
  [`cloud2raster()`](https://georgewoolsey.github.io/cloud2trees/reference/cloud2raster.md)

- Perform individual tree detection using
  [`raster2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/raster2trees.md)

- Quantify individual tree competition metrics using
  [`trees_competition()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_competition.md)
  (*if set to TRUE*)

- Extract tree DBH values from the normalized point cloud using
  [`treels_stem_dbh()`](https://georgewoolsey.github.io/cloud2trees/reference/treels_stem_dbh.md)
  (*if set to TRUE*)

- Model tree DBH values using
  [`trees_dbh()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_dbh.md)
  (*if set to TRUE*)

- Extract tree forest type group using
  [`trees_type()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_type.md)
  (*if set to TRUE*)

- Extract tree CBH values from the normalized point cloud and estimate
  missing values using
  [`trees_cbh()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_cbh.md)
  (*if set to TRUE*)

- Estimate tree biomass (or crown biomass) using
  [`trees_biomass()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_biomass.md)
  (*if method is denoted*)

See the documentation for each individual function called for more
details.

## Usage

``` r
cloud2trees(
  output_dir,
  input_las_dir,
  input_treemap_dir = NULL,
  input_foresttype_dir = NULL,
  input_landfire_dir = NULL,
  accuracy_level = 2,
  max_ctg_pts = 7e+07,
  max_area_m2 = 9e+07,
  transform = FALSE,
  new_crs = NA,
  old_crs = NA,
  keep_intrmdt = FALSE,
  dtm_res_m = 1,
  chm_res_m = 0.25,
  min_height = 2,
  max_height = 70,
  ws = itd_ws_functions()[["log_fn"]],
  estimate_tree_dbh = FALSE,
  max_dbh = 2,
  dbh_model_regional = "cr",
  dbh_model_local = "lin",
  estimate_dbh_from_cloud = FALSE,
  estimate_tree_competition = FALSE,
  competition_buffer_m = 5,
  search_dist_max,
  competition_max_search_dist_m = 10,
  estimate_tree_type = FALSE,
  type_max_search_dist_m = 1000,
  estimate_tree_hmd = FALSE,
  hmd_tree_sample_n = NA,
  hmd_tree_sample_prop = NA,
  hmd_estimate_missing_hmd = FALSE,
  estimate_biomass_method = NA,
  biomass_max_crown_kg_per_m3 = 2,
  estimate_tree_cbh = FALSE,
  cbh_tree_sample_n = NA,
  cbh_tree_sample_prop = NA,
  cbh_which_cbh = "lowest",
  cbh_estimate_missing_cbh = FALSE,
  cbh_min_vhp_n = 3,
  cbh_voxel_grain_size_m = 1,
  cbh_dist_btwn_bins_m = 1,
  cbh_min_fuel_layer_ht_m = 1,
  cbh_lad_pct_gap = 25,
  cbh_lad_pct_base = 25,
  cbh_num_jump_steps = 1,
  cbh_min_lad_pct = 10,
  cbh_frst_layer_min_ht_m = 1,
  overwrite = TRUE
)
```

## Arguments

- output_dir:

  parent directory where new folders `point_cloud_processing_delivery`
  and `point_cloud_processing_temp` will be written for exports

- input_las_dir:

  directory where .las\|.laz point cloud data exists...program will
  search all sub-directories for all .las\|.laz files and process them
  as one

- input_treemap_dir:

  character. directory where Treemap 2016 exists. Use
  [`get_treemap()`](https://georgewoolsey.github.io/cloud2trees/reference/get_treemap.md)
  first.

- input_foresttype_dir:

  character. directory where Forest Type Groups data exists. Use
  [`get_foresttype()`](https://georgewoolsey.github.io/cloud2trees/reference/get_foresttype.md)
  first.

- input_landfire_dir:

  character. directory where LANDFIRE CBD data exists. Use
  [`get_landfire()`](https://georgewoolsey.github.io/cloud2trees/reference/get_landfire.md)
  first.

- accuracy_level:

  numeric. Choose processing accuracy. accuracy_level = 1 uses DTM to
  height normalize the points accuracy_level = 2 uses triangulation with
  high point density (20 pts/m2) to height normalize the points
  accuracy_level = 3 uses triangulation with very high point density
  (100 pts/m2) to height normalize the points

- max_ctg_pts:

  numeric. Max number of points to process at one time. Setting this
  number higher will possibly reduce run times but increase the chance
  of running out of memory and vice versa.

- max_area_m2:

  numeric. Max area to process at one time. See `max_ctg_pts` parameter,
  this one is less important as never experienced memory issues with
  large areas (just lots of points)

- transform:

  logical. should the las/laz files be transformed? If set to `TRUE` the
  parameters `new_crs` must be defined.

- new_crs:

  string. crs to change to as an epsg numerical code

- old_crs:

  string. crs to change from as an epsg numerical code

- keep_intrmdt:

  logical. this process writes intermediate data to the disk, keep those
  intermediate files (classfied, normalized, stem las files)?

- dtm_res_m:

  numeric. The desired resolution of the DTM produced in meters.

- chm_res_m:

  numeric. The desired resolution of the CHM produced in meters.

- min_height:

  numeric. Set the minimum height (m) for individual tree detection

- max_height:

  numeric. Set the maximum height (m) for the canopy height model

- ws:

  numeric or function. Length or diameter of the moving window used to
  detect the local maxima in the units of the input data (usually
  meters). If it is numeric a fixed window size is used. If it is a
  function, the function determines the size of the window at any given
  location on the canopy. By default function takes the height of a
  given pixel as its only argument and return the desired size of the
  search window when centered on that pixel.

- estimate_tree_dbh:

  logical. Should tree DBH be estimated? See
  [`trees_dbh()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_dbh.md).

- max_dbh:

  numeric. Set the largest tree diameter (m) expected in the point cloud

- dbh_model_regional:

  string. Set the model to use for regional dbh-height allometry based
  on FIA tree measurements. Can be "cr" for the Chapman-Richards formula
  (default) or "power" for power function

- dbh_model_local:

  string. Set the model to use for local dbh-height allometry based on
  provided DBH training data in `treels_dbh_locations`. Can be "rf" for
  random forest or "lin" for linear

- estimate_dbh_from_cloud:

  logical. Should DBH be estimated from the point cloud? See
  [`treels_stem_dbh()`](https://georgewoolsey.github.io/cloud2trees/reference/treels_stem_dbh.md).
  Setting to `TRUE` may significantly increase processing time.

- estimate_tree_competition:

  logical. Should tree competition metrics be calculated? See
  [`trees_competition()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_competition.md).
  Setting to `TRUE` may slightly increase processing time.

- competition_buffer_m:

  number. Set buffer around tree (m) to calculate competition metrics

- search_dist_max:

  **\[deprecated\]** Use the `competition_max_search_dist_m` argument
  instead.

- competition_max_search_dist_m:

  number. Maximum search distance (m) to nearest tree for competition.
  Larger search distances will increase processing time and possibly
  result in memory issues. If no competition trees are found within this
  distance, the return column `comp_dist_to_nearest_m` =
  `competition_max_search_dist_m` parameter.

- estimate_tree_type:

  logical. Should tree forest type be estimated? See
  [`trees_type()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_type.md).

- type_max_search_dist_m:

  number. Maximum search distance (m) to obtain forest type group data
  for trees that overlap with non-forest data in the original
  Wilson (2023) data. Larger search distances will increase processing
  time and possibly result in memory issues.

- estimate_tree_hmd:

  logical. Should tree height of the maximum crown diameter (HMD) be
  estimated? See
  [`trees_hmd()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_hmd.md).

- hmd_tree_sample_n, hmd_tree_sample_prop:

  numeric. Provide either `tree_sample_n`, the number of trees, or
  `tree_sample_prop`, the proportion of the trees to attempt to extract
  a HMD from the point cloud for. If neither are supplied,
  `tree_sample_n = 777` will be used. If both are supplied,
  `tree_sample_n` will be used. Increasing `tree_sample_prop` toward
  one (1) will increase the processing time, perhaps significantly
  depending on the number of trees in the `trees_poly` data. The maximum
  number of trees to extract tree HMD using `cloud2trees()` is 20,000.
  Try
  [`trees_hmd()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_hmd.md)
  with outputs from `cloud2trees()` if you want to attempt to extract
  HMD for \>20,000 trees.

- hmd_estimate_missing_hmd:

  logical. It is not likely that HMD will be extracted successfully from
  every tree. Should the missing HMD values be estimated using the tree
  height and location information based on trees for which HMD is
  successfully extracted?

- estimate_biomass_method:

  character. To estimate tree biomass or tree (or crown biomass) enter
  one or a list of multiple biomass methods. See
  [`trees_biomass()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_biomass.md).
  Leave as blank (i.e. `NA`) to skip biomass estimation.

- biomass_max_crown_kg_per_m3:

  numeric. the maximum CBD of the tree crown in kilograms per cubic
  meter. Values above this limit will be set at the median value for the
  area using only stands that have CBD values lower than this limit. The
  default value of 2 kilograms per cubic meter was based on [Mell et al.
  (2009)](https://doi.org/10.1016/j.combustflame.2009.06.015) who found
  the dry bulk density of the tree crown was 2.6 kilograms per cubed
  meter using Douglas-fir trees grown on Christmas tree farms. Set this
  parameter to a large value (e.g. 1e10) or NULL to avoid limiting tree
  crown CBD.

- estimate_tree_cbh:

  logical. Should tree DBH be estimated? See
  [`trees_cbh()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_cbh.md).
  Make sure to set `cbh_estimate_missing_cbh = TRUE` if you want to
  obtain CBH values for cases when CBH cannot be extracted from the
  point cloud.

- cbh_tree_sample_n, cbh_tree_sample_prop:

  numeric. Provide either `tree_sample_n`, the number of trees, or
  `tree_sample_prop`, the proportion of the trees to attempt to extract
  a CBH from the point cloud for. If neither are supplied,
  `tree_sample_n = 333` will be used. If both are supplied,
  `tree_sample_n` will be used. Increasing `tree_sample_prop` toward
  one (1) will increase the processing time, perhaps significantly
  depending on the number of trees in the `trees_poly` data. The maximum
  number of trees to extract tree CBH using `cloud2trees()` is 20,000.
  Try
  [`trees_cbh()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_cbh.md)
  with outputs from `cloud2trees()` if you want to attempt to extract
  CBH for \>20,000 trees.

- cbh_which_cbh:

  character. One of: "lowest"; "highest"; or "max_lad". See Viedma et
  al. (2024) reference.

  - "lowest" - Height of the CBH of the segmented tree based on the last
    distance found in its profile

  - "highest" - Height of the CBH of the segmented tree based on the
    maximum distance found in its profile

  - "max_lad" - Height of the CBH of the segmented tree based on the
    maximum LAD percentage

- cbh_estimate_missing_cbh:

  logical. even if the `cbh_tree_sample_prop` parameter is set to "1",
  it is not likely that CBH will be extracted successfully from every
  tree. Should the missing CBH values be estimated using the tree height
  and location information based on trees for which CBH is successfully
  extracted?

- cbh_min_vhp_n:

  numeric. the minimum number of vertical height profiles (VHPs) needed
  to estimate a CBH.

- cbh_voxel_grain_size_m:

  numeric. horizontal resolution (suggested 1 meter for lad profiles and
  10 meters for LAI maps). See `grain.size` in
  [`leafR::lad.voxels()`](https://rdrr.io/pkg/leafR/man/lad.voxels.html)

- cbh_dist_btwn_bins_m:

  numeric. value for the actual height bin step (in meters). See `step`
  in
  [`LadderFuelsR::get_gaps_fbhs()`](https://rdrr.io/pkg/LadderFuelsR/man/get_gaps_fbhs.html)

- cbh_min_fuel_layer_ht_m:

  numeric. value for the actual minimum base height (in meters). See
  `min_height` in
  [`LadderFuelsR::get_gaps_fbhs()`](https://rdrr.io/pkg/LadderFuelsR/man/get_gaps_fbhs.html)

- cbh_lad_pct_gap:

  numeric. value of the percentile threshold used to identify gaps
  (default percentile 25th). See `perc_gap` in
  [`LadderFuelsR::get_gaps_fbhs()`](https://rdrr.io/pkg/LadderFuelsR/man/get_gaps_fbhs.html)

- cbh_lad_pct_base:

  numeric. value of the percentile threshold used to identify fuels
  layers base height (default percentile 25th). See `perc_base` in
  [`LadderFuelsR::get_gaps_fbhs()`](https://rdrr.io/pkg/LadderFuelsR/man/get_gaps_fbhs.html)

- cbh_num_jump_steps:

  numeric. value for the number of height bin steps that can be jumped
  to reshape fuels layers. See `number_steps` in
  [`LadderFuelsR::get_real_fbh()`](https://rdrr.io/pkg/LadderFuelsR/man/get_real_fbh.html)

- cbh_min_lad_pct:

  numeric. value for the minimum required LAD percentage in a fuel
  layer. See `threshold` in
  [`LadderFuelsR::get_layers_lad()`](https://rdrr.io/pkg/LadderFuelsR/man/get_layers_lad.html)

- cbh_frst_layer_min_ht_m:

  numeric. value for the depth height of the first fuel layer. If the
  first fuel layer has the maximum LAD and its depth is greater than the
  indicated value, then this fuel layer is considered as the CBH of the
  tree. On the contrary, if its depth is \<= the value, the CBH with
  maximum LAD will be the second fuel layer, although it has not the
  maximum LAD. See `hdepth1_height` in
  [`LadderFuelsR::get_cbh_metrics()`](https://rdrr.io/pkg/LadderFuelsR/man/get_cbh_metrics.html)

- overwrite:

  logical. Should the output files in the
  `point_cloud_processing_delivery` directory from previous iterations
  be deleted?

- dbh_model:

  **\[deprecated\]** Use the `dbh_model_regional` or `dbh_model_local`
  argument instead.

## Value

Returns the goods. Exports files of the goods to new folders
"point_cloud_processing_delivery" and "point_cloud_processing_temp" in
the `output_dir` defined by the user in the function call.

## Examples

``` r
 if (FALSE) { # \dontrun{
 # test las file but this could also be a directory path with >1 .las|.laz files
 i <- system.file("extdata", "MixedConifer.laz", package="lidR")
 # run it
 cloud2trees_ans <- cloud2trees::cloud2trees(output_dir = tempdir(), input_las_dir = i)
 # what is it?
 cloud2trees_ans %>% names()
 # there's a DTM
 cloud2trees_ans$dtm_rast %>% terra::plot()
 # there's a CHM
 cloud2trees_ans$chm_rast %>% terra::plot()
 # there are tree crowns
 cloud2trees_ans$crowns_sf %>% dplyr::glimpse()
 cloud2trees_ans$crowns_sf %>% ggplot2::ggplot() +
  ggplot2::geom_sf(mapping = ggplot2::aes(fill = tree_height_m))
 # there are tree top points
 cloud2trees_ans$treetops_sf %>% dplyr::glimpse()
 cloud2trees_ans$treetops_sf %>% ggplot2::ggplot() +
  ggplot2::geom_sf(mapping = ggplot2::aes(color = tree_height_m))
 } # }
```
