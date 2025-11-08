# Download url data

Generic function to download data from a url with .zip file data
Functions
[`get_treemap()`](https://georgewoolsey.github.io/cloud2trees/reference/get_treemap.md)
and
[`get_foresttype()`](https://georgewoolsey.github.io/cloud2trees/reference/get_foresttype.md)
use this

## Usage

``` r
get_url_data(
  eval_url = NULL,
  my_name = NULL,
  savedir = NULL,
  req_file_list = NULL,
  force = F,
  cleanup_zip = T,
  move_files_to_top = T
)
```

## Arguments

- eval_url:

  Required url of the .zip file to be downloaded. url string must end in
  .zip

- my_name:

  Required name for the folder which data will be extracted to

- savedir:

  Optional directory to save data in a new location. Defaults to package
  contents.

- req_file_list:

  Optional list of files to check for before re-downloading full data

- force:

  Whether to overwrite existing data

- cleanup_zip:

  Whether to remove the .zip file after extracting the contents

## Examples

``` r
 if (FALSE) { # \dontrun{
 get_treemap()
 } # }
```
