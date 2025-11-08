# Tile raw `.las`\|`.laz` files to work with smaller chunks

Function to tile raw `.las`\|`.laz` files to work with smaller chunks
based on point density and coverage area

## Usage

``` r
chunk_las_catalog(
  folder,
  outfolder = getwd(),
  accuracy_level = 2,
  max_ctg_pts = 7e+07,
  max_area_m2 = 9e+07,
  transform = FALSE,
  new_crs = NA,
  old_crs = NA
)
```

## Arguments

- folder:

  string. The path of a folder containing a set of las/laz files. Can
  also be a vector of file paths.

- outfolder:

  string. The path of a folder to write the tiled las files to.

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

## Value

A list of 1) `process_data` an sf object; 2) `is_chunked_grid` indicator
if chunks were created; 3) `plt` a ggplot object

## References

<https://r-lidar.github.io/lidRbook/norm.html>
<https://github.com/r-lidar/lasR/issues/18#issuecomment-2027818414>

## Examples

``` r
 if (FALSE) { # \dontrun{
 f <- "../lasdata"
 chunk_las_catalog(folder = f, outfolder = getwd())
 } # }
```
