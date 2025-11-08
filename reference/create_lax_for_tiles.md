# Create spatial index `.lax` files

Function to create spatial index files .lax for .las\|.laz files to
speed up processing

## Usage

``` r
create_lax_for_tiles(las_file_list)
```

## Arguments

- las_file_list:

  a list of .las\|.laz files with full directory path

## Value

A list of file names

## References

<https://r-lidar.github.io/lidRbook/spatial-indexing.html>

## Examples

``` r
 if (FALSE) { # \dontrun{
 f <- list.files(getwd(), pattern = ".*\\.(laz|las)$", full.names = TRUE)
 create_lax_for_tiles(las_file_list = f)
 } # }
```
