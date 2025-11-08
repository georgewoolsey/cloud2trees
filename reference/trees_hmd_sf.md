# Estimate HMD using tree crown polygons and normalized point cloud data

`trees_hmd_sf()` uses the input tree crown polygons (e.g. as exported by
[`raster2trees()`](https://georgewoolsey.github.io/cloud2trees/reference/raster2trees.md))
with the columns `treeID` and `tree_height_m` to extracting the height
of the maximum crown diameter (HMD) using height normalized point cloud
data (e.g. as exported by
[`cloud2raster()`](https://georgewoolsey.github.io/cloud2trees/reference/cloud2raster.md)).

HMD is extracted directly from the height normalized point cloud by
finding the height of the non-ground point farthest from the tree center
(i.e. tree top).

An early version of this process was developed by [Andrew Sanchez
Meador](https://github.com/bi0m3trics).

There are likely to be trees for which there is insufficient data in the
point cloud to successfully estimate HMD. The user can elect to estimate
missing HMD values which is accomplished via:

- Attempt to extract HMD from all trees

- Successfully extracted HMD trees become training data used to estimate
  the height-HMD allometry relationship that is spatially informed using
  the relative tree location compared to the training data

- The height and location predicting HMD model built from the point
  cloud training data is used to predict HMD for the non-training (i.e.
  missing HMD) data

## Usage

``` r
trees_hmd_sf(
  trees_poly,
  norm_las = NULL,
  tree_sample_n = NA,
  tree_sample_prop = NA,
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
  .laz\|.las file", -or- an object of class `LAScatalog`. It is your
  responsibility to ensure that the point cloud is projected the same as
  the `trees_poly` data

- tree_sample_n, tree_sample_prop:

  numeric. Provide either `tree_sample_n`, the number of trees, or
  `tree_sample_prop`, the proportion of the trees to attempt to extract
  a HMD from the point cloud for. If neither are supplied,
  `tree_sample_n = 777` will be used. If both are supplied,
  `tree_sample_n` will be used. Increasing `tree_sample_prop` toward
  one (1) will increase the processing time, perhaps significantly
  depending on the number of trees in the `trees_poly` data.

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
  prefixed with "hmd\_". Leave as NA to return a data.frame of the trees
  from tree list from `trees_poly` with HMD values added
