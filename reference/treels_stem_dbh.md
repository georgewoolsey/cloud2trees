# detect tree stems and estimate DBH using `TreeLS` package

`treels_stem_dbh()` is an all-in-one function to process height
normalized .las\|.laz files using the functionality of the `TreeLS`
package. The function generates a list of stems and estimates DBH
directly from the point cloud. `treels_stem_dbh()` outputs:

- A .laz file with the `Classification` data updated to: ground points
  (class 2); water points (class 9); stem points (class 4); non-stem
  (class 5).

- A vector data file in `gpkg` format with the tree identification stem
  locations, heights, and DBH estimates.

The order of operations is:

- Detect tree stems/boles from the height normalized point cloud using
  [`TreeLS::treeMap()`](https://rdrr.io/pkg/TreeLS/man/treeMap.html)
  with the
  [`TreeLS::map.hough()`](https://rdrr.io/pkg/TreeLS/man/map.hough.html)
  algorithm

- Merge overlapping tree coordinates using
  [`TreeLS::treeMap.merge()`](https://rdrr.io/pkg/TreeLS/man/treeMap.merge.html)

- Assign tree IDs to the original points using
  [`TreeLS::treePoints()`](https://rdrr.io/pkg/TreeLS/man/treePoints.html)
  with the
  [`TreeLS::trp.crop()`](https://rdrr.io/pkg/TreeLS/man/trp.crop.html)
  algorithm

- Flag only the stem points using
  [`TreeLS::stemPoints()`](https://rdrr.io/pkg/TreeLS/man/stemPoints.html)
  with the
  [`TreeLS::stm.hough()`](https://rdrr.io/pkg/TreeLS/man/stm.hough.html)
  algorithm

- Perform DBH estimation using
  [`TreeLS::tlsInventory()`](https://rdrr.io/pkg/TreeLS/man/tlsInventory.html)
  with the
  [`TreeLS::shapeFit()`](https://rdrr.io/pkg/TreeLS/man/shapeFit.html)
  algorithm

## Usage

``` r
treels_stem_dbh(
  folder,
  outfolder,
  min_height = 2,
  max_dbh = 2,
  chunk_these = FALSE
)
```

## Arguments

- folder:

  string. The path of a folder containing a set of las/laz files. Can
  also be a vector of file paths.

- outfolder:

  string. The path of a folder to write the tiled vector files to

- min_height:

  numeric. Set the minimum height (m) for individual tree detection

- max_dbh:

  numeric. Set the largest tree diameter (m) expected in the point cloud

- chunk_these:

  logical. Do the las/laz files need to be tiled to work with smaller
  subsets? See `is_chunked_grid` in
  [`chunk_las_catalog()`](https://georgewoolsey.github.io/cloud2trees/reference/chunk_las_catalog.md)

## Value

Returns an `sf` data.frame with TreeLS detected trees and DBH estimated
directly from the point cloud. Exports files in the `outfolder` defined
by the user in the function call.

## References

<https://github.com/tiagodc/TreeLS>

## Examples

``` r
 if (FALSE) { # \dontrun{
 o <- "../data"
 i <- "../data/normlasdata"
 r <- cloud2trees::treels_stem_dbh(folder = i, outfolder = o)
 r %>% names()
 } # }
```
