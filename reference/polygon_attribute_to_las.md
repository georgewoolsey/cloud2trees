# function to attach polygon attribute to point cloud

`polygon_attribute_to_las()` function to attach polygon attribute to
point cloud

## Usage

``` r
polygon_attribute_to_las(las, poly_df, attribute, force_crs = F)
```

## Arguments

- las:

  an object of class LAS

- poly_df:

  an object of class sf with only POLYGON geometry (use
  cloud2trees::simplify_multipolygon_crowns() first)

- attribute:

  character. a data attribute in the `poly_df` that you want to
  spatially attach to the `las` (e.g. treeID)

- force_crs:

  logical. turn on the force same crs parameter if confident that data
  are in same projection.

## Value

Returns an `LAS` class object

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
 } # }
```
