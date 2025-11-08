# Simplify MULTIPOLYGON to POLYGON geometry in an `sf` class object

Function to simplify MULTIPOLYGON geometry to POLYGON geometry in an
`sf` class object by selecting the largest segment of the MULTIPOLYGON

## Usage

``` r
simplify_multipolygon_crowns(trees_poly)
```

## Arguments

- trees_poly:

  data.frame. A data frame of `sf` class with POLYGON,MULTIPOLYGON
  geometry (see
  [`sf::st_geometry_type()`](https://r-spatial.github.io/sf/reference/st_geometry_type.html))
  and the column `treeID`

## Value

A `sf` class object data frame

## Examples

``` r
 if (FALSE) { # \dontrun{
 f <- paste0(system.file(package = "cloud2trees"),"/extdata/crowns_poly.gpkg")
 crowns <- sf::st_read(f, quiet = T)
 crowns %>% sf::st_geometry_type() %>% table()
 crowns_simp <- simplify_multipolygon_crowns(crowns)
 crowns_simp %>% sf::st_geometry_type() %>% table()
 } # }
```
