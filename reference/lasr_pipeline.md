# Create `lasR` package pipeline to process a las grid tile created via [`chunk_las_catalog()`](https://georgewoolsey.github.io/cloud2trees/reference/chunk_las_catalog.md)

Create `lasR` package pipeline to process a las grid tile created via
[`chunk_las_catalog()`](https://georgewoolsey.github.io/cloud2trees/reference/chunk_las_catalog.md)

## Usage

``` r
lasr_pipeline(
  processing_grid_num = 1,
  process_data,
  keep_intrmdt = F,
  dtm_res_m = 1,
  chm_res_m = 0.25,
  min_height = 2,
  max_height = 70,
  dtm_dir = getwd(),
  chm_dir = getwd(),
  classify_dir = getwd(),
  normalize_dir = getwd()
)
```

## Arguments

- processing_grid_num:

  numeric. processing_grid column in the data.frame created via
  [`chunk_las_catalog()`](https://georgewoolsey.github.io/cloud2trees/reference/chunk_las_catalog.md)

- process_data:

  data.frame. data.frame created via
  [`chunk_las_catalog()`](https://georgewoolsey.github.io/cloud2trees/reference/chunk_las_catalog.md)

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

- dtm_dir:

  string. The path of a folder to write the tiled DTM files to.

- chm_dir:

  string. The path of a folder to write the tiled CHM files to.

- classify_dir:

  string. The path of a folder to write the classified .las files to.

- normalize_dir:

  string. The path of a folder to write the normalized .las files to.

## Value

A `lasR` pipeline answer list

## References

<https://r-lidar.github.io/lasR/index.html>
