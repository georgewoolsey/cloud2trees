# Download Forest Type Groups of the Continental United States data

The LANDFIRE Forest Canopy Bulk Density (CBD) data is used to estimate
individual tree crown biomass in kilograms See
[`trees_biomass_landfire()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_biomass_landfire.md)

## Usage

``` r
get_landfire(savedir = NULL, force = F)
```

## Arguments

- savedir:

  Optional directory to save data in a new location. Defaults to package
  contents.

- force:

  Whether to overwrite existing data

## References

- [LANDFIRE Forest Canopy Bulk Density
  (CBD)](https://landfire.gov/fuel/cbd) U.S. Department of Agriculture
  and U.S. Department of the Interior.

## Examples

``` r
 if (FALSE) { # \dontrun{
 get_landfire()
 } # }
```
