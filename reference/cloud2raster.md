# Use raw .las\|.laz files to generate CHM, DTM, and Normalized .las files

`cloud2raster()` is an all-in-one function to process raw .las\|.laz
files to generate a CHM raster (.tif), a DTM raster (.tif), and .las
files which have been height normalized. The order of operations is:

- Tile the raw point cloud to work with smaller chunks and reduce the
  potential for memory issues with high density clouds using
  [`chunk_las_catalog()`](https://georgewoolsey.github.io/cloud2trees/reference/chunk_las_catalog.md)

- Classify the point cloud using
  [`lasR::classify_with_csf()`](https://rdrr.io/pkg/lasR/man/classify_with_csf.html)

- Remove outlier points using
  [`lasR::classify_with_ivf()`](https://rdrr.io/pkg/lasR/man/classify_with_ivf.html)

- Produce a triangulation of the ground points (meshed DTM) using
  [`lasR::triangulate()`](https://rdrr.io/pkg/lasR/man/triangulate.html)

- Rasterize the result of the Delaunay triangulation using
  [`lasR::rasterize()`](https://rdrr.io/pkg/lasR/man/rasterize.html) to
  create a DTM

- Height normalize the point cloud using either the DTM or the
  triangulation
  [`lasR::transform_with()`](https://rdrr.io/pkg/lasR/man/transform_with.html)

- Use the height normalized point cloud to create the CHM based on the
  highest point in a pixel using
  [`lasR::rasterize()`](https://rdrr.io/pkg/lasR/man/rasterize.html)

- Pits and spikes filling of the CHM raster using
  [`lasR::pit_fill()`](https://rdrr.io/pkg/lasR/man/pit_fill.html)

- Smooth the CHM raster tile gaps using
  [`terra::focal()`](https://rspatial.github.io/terra/reference/focal.html)

## Usage

``` r
cloud2raster(
  output_dir,
  input_las_dir,
  input_treemap_dir = NULL,
  input_foresttype_dir = NULL,
  accuracy_level = 2,
  max_ctg_pts = 7e+07,
  max_area_m2 = 9e+07,
  transform = FALSE,
  new_crs = NA,
  old_crs = NA,
  keep_intrmdt = F,
  dtm_res_m = 1,
  chm_res_m = 0.25,
  min_height = 2,
  max_height = 70,
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

- overwrite:

  logical. Should the output files in the
  `point_cloud_processing_delivery` directory from previous iterations
  be deleted?

## Value

Returns the goods. Exports files of the goods to new folders
"point_cloud_processing_delivery" and "point_cloud_processing_temp" in
the `output_dir` defined by the user in the function call.

## References

<https://r-lidar.github.io/lasR/index.html>
<https://r-lidar.github.io/lidRbook/normalization.html>

## Examples

``` r
 if (FALSE) { # \dontrun{
 # test las file but this could also be a directory path with >1 .las|.laz files
 i <- system.file("extdata", "MixedConifer.laz", package="lidR")
 # run it
 r <- cloud2trees::cloud2raster(output_dir = tempdir(), input_las_dir = i)
 # what is it?
 r %>% names()
 # there's a DTM
 r$dtm_rast %>% terra::plot()
 # there's a CHM
 r$chm_rast %>% terra::plot()
 # there's a data.frame with the file structure for the project
 r$create_project_structure_ans %>% dplyr::glimpse()
 # there's a information detailing how the point cloud was processed
 r$chunk_las_catalog_ans$process_data %>% dplyr::glimpse()
 r$chunk_las_catalog_ans$is_chunked_grid
 r$chunk_las_catalog_ans$las_ctg@data %>% dplyr::glimpse()
 # there's a list of the height normalized .las files created
 r$normalize_flist
 } # }
```
