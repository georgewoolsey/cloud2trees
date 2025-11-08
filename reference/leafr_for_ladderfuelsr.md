# re-writes `leafR` steps to allow for treeID as input for `ladderfuelsR` or [`ladderfuelsr_cbh()`](https://georgewoolsey.github.io/cloud2trees/reference/ladderfuelsr_cbh.md)

`leafr_for_ladderfuelsr()` is a re-write of
[`leafR::lad.voxels()`](https://rdrr.io/pkg/leafR/man/lad.voxels.html)
and
[`leafR::lad.profile()`](https://rdrr.io/pkg/leafR/man/lad.profile.html)
to:

- removes the requirement to use a file written to disk

- allows for the calculation of LAD by the `treeID` attribute so that
  don't have to pass individual tree point clouds

- updates to the use of the latest `lidR` functionality and removes the
  use of `sp` and `raster` functions

- updates the function to `tidy` data manipulation

## Usage

``` r
leafr_for_ladderfuelsr(
  las,
  voxel_grain_size_m = 1,
  k = 1,
  attribute = "treeID",
  min_pulses = 0,
  relative = FALSE
)
```

## Arguments

- las:

  an object of class LAS that has been height normalized.

- voxel_grain_size_m:

  numeric. horizontal resolution (suggested 1 meter for lad profiles and
  10 meters for LAI maps). See `grain.size` in
  [`leafR::lad.voxels()`](https://rdrr.io/pkg/leafR/man/lad.voxels.html)

- k:

  numeric. coefficient to transform effective LAI to real LAI (k = 1;
  for effective LAI)

- attribute:

  character. The column name of the attribute to group the return LAD
  profile by. Default is "treeID". The attribute (whatever it is defined
  as) must exist in the `las` data

- min_pulses:

  numeric. minimum number of pulses required to return a record by
  attribute. set to zero (default) to leave data unfiltered.

- relative:

  logical. produce lad profile by relative total LAI values. Indicate
  when using effective LAI value. if set to TRUE, lad value will be
  relative_lad; otherwise, lad value will be mean_lad

## Value

Returns an data.frame which can have multiple treeIDs to use as input
for
[`ladderfuelsr_cbh()`](https://georgewoolsey.github.io/cloud2trees/reference/ladderfuelsr_cbh.md)

## References

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
 } # }
```
