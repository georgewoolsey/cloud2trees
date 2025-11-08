# Check a data frame for spatial point data. Convert to points if needed.

Check a data frame for spatial point data. Convert to points if needed.

## Usage

``` r
check_spatial_points(tree_list, crs = NA)
```

## Arguments

- tree_list:

  data.frame. A data frame with the columns `treeID`, `tree_x`,
  `tree_y`, and `tree_height_m`. If an `sf` class object with POINT
  geometry (see
  [`sf::st_geometry_type()`](https://r-spatial.github.io/sf/reference/st_geometry_type.html)),
  the program will use the data "as-is" and only require the `treeID`
  and `tree_height_m` columns.

- crs:

  string. A crs string as returned from
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)
  or the EPSG code of the x,y coordinates. Defaults to the crs of the
  `tree_list` data if of class "sf".
