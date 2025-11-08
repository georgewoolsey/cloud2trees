# Find the location of external data

Find the location of external data Functions `get_*()` download external
data

## Usage

``` r
find_ext_data(
  input_treemap_dir = NULL,
  input_foresttype_dir = NULL,
  input_landfire_dir = NULL
)
```

## Arguments

- input_treemap_dir:

  character. directory where Treemap 2016 exists. Use
  [`get_treemap()`](https://georgewoolsey.github.io/cloud2trees/reference/get_treemap.md)
  first.

- input_foresttype_dir:

  character. directory where Forest Type Groups data exists. Use
  [`get_foresttype()`](https://georgewoolsey.github.io/cloud2trees/reference/get_foresttype.md)
  first.

- input_landfire_dir:

  character. directory where LANDFIRE CBD data exists. Use
  [`get_landfire()`](https://georgewoolsey.github.io/cloud2trees/reference/get_landfire.md)
  first.

## Value

Returns a list where the values will be either NULL if unable to locate
the external data files , or the directory where the external data files
were located. The list includes the named variables `treemap_dir` and
`foresttype_dir`

## Examples

``` r
 if (FALSE) { # \dontrun{
 find_ext_data()
 } # }
```
