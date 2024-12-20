# cloud2trees 0.3.1

Updates the processes that rely on the `lasR` package which implemented major revisions with breaking changes in [`lasR 0.13.0`](https://r-lidar.github.io/lasR/news/index.html#lasr-0130). The package `cloud2trees` now requires `lasR >= 0.13.1`. Even with these updates, users may continue to encounter some bugs in the near future. See (@georgewoolsey, #3).

To update execute:

```r
# update lasR
install.packages("lasR", repos = "https://r-lidar.r-universe.dev")
# update cloud2trees
remotes::install_github(repo = "georgewoolsey/cloud2trees", upgrade = F)
```

- Change: `lasr_*()` (internal functions) replace `LASlib` filters to conform with new `lasR` filtering (see [https://r-lidar.github.io/lasR/reference/filters.html]())
- Change: `lasr_dtm_norm()` (internal function) replaces the `LASlib` filter `-keep_random_fraction` with `lasR::sampling_pixel()` to decimate the ground points (see [https://github.com/r-lidar/lasR/issues/102]())

# cloud2trees 0.3.0

Implements a process to extract the USDA Forest Inventory and Analysis (FIA) forest type group based on the [Forest Type Groups of the Continental United States](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d) data (Wilson 2023). See (@georgewoolsey, #2).

- New: Adds the function `trees_type()` to use the input tree list (e.g. as exported by `raster2trees()`) to attach species information using FIA codes
- New: Adds `get_data()` as all-in-one function that downloads all of the external data used by the package
- New: Adds `get_foresttype()` downloads the external data [Forest Type Groups of the Continental United States](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d)
- New: Adds `find_ext_data()` attempts to find the location of external data (as downloaded from `get_*()`)
- Change: `cloud2trees()` incorporates the CONUS forest type process via `trees_type()`

# cloud2trees 0.2.0

Integrates a process to extract the crown base height (CBH) using the workflow outlined in [Viedma et al. (2024)](https://doi.org/10.1111/2041-210X.14427) and using their package `LadderFuelsR` ([https://github.com/olgaviedma/LadderFuelsR]()). See (@georgewoolsey, #1).

- New: Adds the function `trees_cbh()` to use input tree crown polygons (e.g. as exported by `raster2trees()`) to estimate tree CBH based on an input normalized point cloud
- New: Adds `ladderfuelsr_cbh()` as all-in-one function to extract CBH from a single tree, height normalized cloud using the functionality of the `LadderFuelsR` package
- Change: `cloud2trees()` incorporates the CBH process via `trees_cbh()`
- New: Adds error handling in `cloud2trees()` for stages after trees are extracted (e.g. during DBH extraction) the function will complete successfully and write the extracted trees data to the `point_cloud_processing_delivery` directory. An error message will be issued at the end (*ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! in: ...*) and state which portion did not complete successfully so that the user can attempt to re-run that portion with different settings (generally a `trees_*()` function)

# cloud2trees 0.1.0

- Open to public
