# estimate CBH for a *single tree* using `LadderFuelsR` package

`ladderfuelsr_cbh()` is an all-in-one function to process height
normalized .las\|.laz files using the functionality of the
`LadderFuelsR` package. The function returns a a list of data.frame
objects that are the results of the different LadderFuelsR steps.
Returns NULL if the process is unable to detect a CBH from the point
cloud. `ladderfuelsr_cbh()` outputs:

The order of operations is:

- Create a data frame of the 3D voxels information (xyz) with Leaf Area
  Density (LAD) values from las file using
  [`leafR::lad.voxels()`](https://rdrr.io/pkg/leafR/man/lad.voxels.html)

- Calculate the lad profile from the input lad.voxels using
  [`leafR::lad.profile()`](https://rdrr.io/pkg/leafR/man/lad.profile.html)

- Calculate gaps and fuel layers base height (FBH) as the difference in
  percentiles between consecutive LAD values along the vertical tree
  profile (VTP) using
  [`LadderFuelsR::get_gaps_fbhs()`](https://rdrr.io/pkg/LadderFuelsR/man/get_gaps_fbhs.html)

- Calculate the percentile value of each fuel layer base height using
  [`LadderFuelsR::calculate_gaps_perc()`](https://rdrr.io/pkg/LadderFuelsR/man/calculate_gaps_perc.html)

- Calculate distances (and their heights) between fuel layers as the
  difference between consecutive gaps and fuel bases using
  [`LadderFuelsR::get_distance()`](https://rdrr.io/pkg/LadderFuelsR/man/get_distance.html)

- Calculate fuels depth as the difference between gaps interleaved
  between fuel layers minus one step if the fuel depths are greater than
  one step using
  [`LadderFuelsR::get_depths()`](https://rdrr.io/pkg/LadderFuelsR/man/get_depths.html)

- Reshape fuel layers after removing distances equal to any number of
  height bin steps, keeping the first "base height" from those
  consecutive ones separated by such distance using
  [`LadderFuelsR::get_real_fbh()`](https://rdrr.io/pkg/LadderFuelsR/man/get_real_fbh.html)

- Recalculate fuel layers depth after considering distances greater than
  the actual height bin step using
  [`LadderFuelsR::get_real_depths()`](https://rdrr.io/pkg/LadderFuelsR/man/get_real_depths.html)

- Recalculate the distance between fuel layers after considering
  distances greater than any number of height bin steps using
  [`LadderFuelsR::get_effective_gap()`](https://rdrr.io/pkg/LadderFuelsR/man/get_effective_gap.html)

- Calculate the percentage of LAD within each fuel layer (first output)
  and removes those fuel layers with LAD percentage less than a
  specified threshold using
  [`LadderFuelsR::get_layers_lad()`](https://rdrr.io/pkg/LadderFuelsR/man/get_layers_lad.html)

- Determine the CBH of a segmented tree using three criteria: maximum
  LAD percentage, maximum distance and the last distance using
  [`LadderFuelsR::get_cbh_metrics()`](https://rdrr.io/pkg/LadderFuelsR/man/get_cbh_metrics.html)

## Usage

``` r
ladderfuelsr_cbh(
  lad_profile_df = NULL,
  las = NULL,
  treeID = NA,
  min_vhp_n = 4,
  voxel_grain_size_m = 1,
  dist_btwn_bins_m = 1,
  min_fuel_layer_ht_m = 1,
  lad_pct_gap = 25,
  lad_pct_base = 25,
  num_jump_steps = 1,
  min_lad_pct = 10,
  frst_layer_min_ht_m = 1
)
```

## Arguments

- lad_profile_df:

  data.frame. the return of
  [`leafr_for_ladderfuelsr()`](https://georgewoolsey.github.io/cloud2trees/reference/leafr_for_ladderfuelsr.md)
  or a data.frame that must have the columns: treeID, lad, total_pulses,
  height. if both the `lad_profile_df` and `las` parameters are defined,
  the preference is the data.frame

- las:

  string -or- object. a single tree .las\|.laz file path -OR- an object
  of class LAS that has been height normalized.

- treeID:

  numeric. the LadderFuelsR process requires a treeID that uniquely
  identifies points within a tree , if left as NA this process will
  attempt to locate the `treeID` data based on an attribute in the point
  cloud or data.frame which should be numeric

- min_vhp_n:

  numeric. the minimum number of vertical height profiles (VHPs) needed
  to estimate a CBH.

- voxel_grain_size_m:

  numeric. only used if `las` parameter is defined. horizontal
  resolution (suggested 1 meter for lad profiles). See `grain.size` in
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
  indicated value , then this fuel layer is considered as the CBH of the
  tree. On the contrary, if its depth is \<= the value, the CBH with
  maximum LAD will be the second fuel layer, although it has not the
  maximum LAD. See `hdepth1_height` in
  [`LadderFuelsR::get_cbh_metrics()`](https://rdrr.io/pkg/LadderFuelsR/man/get_cbh_metrics.html)

## Value

Returns an list of data.frame objects that are the results of the
different LadderFuelsR steps. Returns NULL if the process is unable to
detect a CBH from the point cloud.

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
 # polygon data
 f <- system.file(package = "cloud2trees","extdata","crowns_poly.gpkg")
 trees_poly <- sf::st_read(f)
 # simplify polygons
 trees_poly <- simplify_multipolygon_crowns(trees_poly)
 # point cloud data
 lf <- system.file(package = "cloud2trees","extdata","norm_las","RMNP_017_2018_normalize.las")
 las <- lidR::readLAS(lf)
 las@data %>% dplyr::glimpse()
 # polygon_attribute_to_las to attach treeID to las
 las <- polygon_attribute_to_las(las, trees_poly, force_crs = T, attribute = "treeID")
 las@data %>% dplyr::glimpse()
 # get the lad profile for each treeID
 lad_profile <- leafr_for_ladderfuelsr(
     las
     , voxel_grain_size_m = 1
     , k = 1
     , group_treeID = T
     , relative = F
   )
 dplyr::glimpse(lad_profile)
 # extract the CBH using ladderfuelsr_cbh()
 # before we extract the CBH using ladderfuelsr_cbh(), treeID has to be numeric
 lad_profile <- lad_profile %>%
   dplyr::mutate(
     treeID_backup = treeID, treeID = as.factor(treeID)
   )
 # for one tree
 ladderfuelsr_cbh(
     lad_profile_df = lad_profile
     , treeID = lad_profile$treeID[1]
   ) %>%
   purrr::pluck("cbh_metrics") %>%
   dplyr::glimpse()
 # we can map over multiple trees
 cbhs <- lad_profile$treeID %>%
   unique() %>%
   .[1:22] %>%
   purrr::map(\(x) ladderfuelsr_cbh(
       lad_profile_df = lad_profile
       , treeID = x
     ) %>%
     purrr::pluck("cbh_metrics")
   ) %>%
   dplyr::bind_rows()
 dplyr::glimpse(cbhs)
 ggplot2::ggplot(data = cbhs, mapping = ggplot2::aes(x=last_Hcbh)) +
   ggplot2::geom_density()
 } # }
```
