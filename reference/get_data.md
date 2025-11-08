# Download all external data used by package

The package requires external data to estimate individual tree DBH (see
[`get_treemap()`](https://georgewoolsey.github.io/cloud2trees/reference/get_treemap.md))
and to extract the forest type for a tree list (see
[`get_foresttype()`](https://georgewoolsey.github.io/cloud2trees/reference/get_foresttype.md)).
This is a all-in-one function that downloads all of the external data
used by the package.

## Usage

``` r
get_data(savedir = NULL, force = F)
```

## Arguments

- savedir:

  Optional directory to save data in a new location. Defaults to package
  contents.

- force:

  Whether to overwrite existing data

## Examples

``` r
 if (FALSE) { # \dontrun{
 get_data()
 } # }
```
