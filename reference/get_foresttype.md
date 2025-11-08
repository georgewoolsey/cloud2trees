# Download Forest Type Groups of the Continental United States data

The Forest Type Groups of the Continental United States data is used to
estimate individual tree forest type group. See
[`trees_type()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_type.md)

## Usage

``` r
get_foresttype(savedir = NULL, force = F, res = 30)
```

## Arguments

- savedir:

  Optional directory to save data in a new location. Defaults to package
  contents.

- force:

  Whether to overwrite existing data

- res:

  Resolution of forest type data to download. Default of "30" downloads
  30m raster, any other value will download 90m raster.

## References

- [Forest Type Groups of the Continental United
  States](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d)
  Wilson, B.T. (2023). Forest Type Groups of the Continental United
  States.

## Examples

``` r
 if (FALSE) { # \dontrun{
 get_foresttype()
 } # }
```
