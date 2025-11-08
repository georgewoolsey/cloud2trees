# Create project structure

Function to generate nested project directories based on the
user-defined directory to create output file structure

## Usage

``` r
create_project_structure(
  output_dir,
  input_las_dir,
  input_treemap_dir = NULL,
  input_foresttype_dir = NULL
)
```

## Arguments

- output_dir:

  parent directory where new folders `point_cloud_processing_delivery`
  and `point_cloud_processing_temp` will be written for exports

- input_las_dir:

  directory where .las\|.laz point cloud data exists...program will
  search all sub-directories for all .las\|.laz files and process them
  as one

- input_treemap_dir:

  character. directory where Treemap 2016 exists. Use
  [`get_treemap()`](https://georgewoolsey.github.io/cloud2trees/reference/get_treemap.md)
  first.

- input_foresttype_dir:

  character. directory where Forest Type Groups data exists. Use
  [`get_foresttype()`](https://georgewoolsey.github.io/cloud2trees/reference/get_foresttype.md)
  first.

## Value

A data.frame.
