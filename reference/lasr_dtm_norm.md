# use lasR to combine DTM and normalize step

Combining DTM and normalize step using `lasR` functionality

## Usage

``` r
lasr_dtm_norm(dtm_file_name, frac_for_tri = 1, dtm_res = 1, norm_accuracy = 2)
```

## Arguments

- dtm_file_name:

  string. Where to write the DTM.

- frac_for_tri:

  numeric. The fraction of points used in Delauny triangulation.

- dtm_res:

  numeric. The desired resolution of the DTM produced in meters.

- norm_accuracy:

  numeric. see `chunk_las_catalog`. Choose processing accuracy.
  accuracy_level = 1 uses DTM to height normalize the points
  accuracy_level = 2 uses triangulation with high point density (20
  pts/m2) to height normalize the points accuracy_level = 3 uses
  triangulation with very high point density (100 pts/m2) to height
  normalize the points

## Value

A `lasR` pipeline

## References

<https://r-lidar.github.io/lidRbook/normalization.html>
<https://github.com/r-lidar/lasR/issues/18#issuecomment-2027818414>
