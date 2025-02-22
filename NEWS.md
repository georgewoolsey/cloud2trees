# cloud2trees 0.5.7

# cloud2trees 0.5.6

- Fix: `raster2trees()` would potentially fail when writing data if `tempdir = tempdir()` and the raster was too big to fit in memory. This update also forces the `treeID` column in the return data to character type. Lastly, output from the function is now written as "final_detected_\*.gpkg" instead of "chm_detected_\*.gpkg" to match the output from `cloud2trees()`.
- Change: `cloud2trees()` now writes the return data from `raster2trees()` ("final_detected_\*.gpkg") prior to completing the following steps and the overwrites this file upon completion of the other steps.

# cloud2trees 0.5.5

- Fix: `raster2trees()` would potentially fail when processing large rasters that where not read by `terra` as "in-memory" due to invalid tree crown geometries generated during the raster tile processing. This update implements additional checks and fixes in the processing section when the raster is too big to fit in memory.

# cloud2trees 0.5.4

The `trees_biomass()` (and `trees_biomass_*()`) function allowed for the application of unconstrained tree crown bulk density (CBD) values in kilograms per cubic meter to calculate crown biomass in kilograms. These CBD values are calculated using the sequence of equations detailed in the `trees_biomass_*()` function documentation. In scenarios where there were only a few small trees with small crown diameters and short crown lengths, for example, CBD values of >5 kilograms per cubic meter (sometimes even much larger) were estimated using this process. These high CBD values are improbable based on the literature. [Mell et al. (2009)](https://doi.org/10.1016/j.combustflame.2009.06.015) found the dry bulk density of the tree crown was 2.6 kilograms per cubed meter using Douglas-fir trees grown on Christmas tree farms. This update allows for users to constrain tree CBD values.

- Change: `trees_biomass()` (and the `trees_biomass_*()` functions) now have a parameter `max_crown_kg_per_m3` which limits the maximum CBD of the tree crown in kilograms per cubic meter. Values above this limit will be set at the median value for the area using only stands that have CBD values lower than this limit.
- Change: `cloud2trees()` adds the `biomass_max_crown_kg_per_m3` parameter to limit the maximum CBD of the tree crown in kilograms per cubic meter.

# cloud2trees 0.5.3

- New: Adds the function `itd_tuning()` to to visually assess tree crown delineation results from different window size functions used for the detection of individual trees.
- Change: `cloud2trees()` and `raster2trees()` now use the same non-linear window function for individual tree detection (at some point the settings got out of alignment)

# cloud2trees 0.5.2

- Fix: `trees_cbh()` could not check CBH against height and could not estimate missing CBH values if an attribute named "treeID" pre-existed in the point cloud data. This fix forces the overwrite of the "treeID" attribute in the point cloud data using polygon data.
- Change: `cloud2trees()` limits the number of trees to extract CBH from to 20,000. If users desire to attempt to extract CBH for >20,000 trees, `trees_cbh()` can be used in standalone with outputs from `cloud2trees()`.

# cloud2trees 0.5.1

Updates the process to extract tree CBH by implementing a re-creation of some of the steps needed that were initially developed in other packages. Improves performance of extraction of CBH from the point cloud to possibly allow users to increase the `tree_sample_n` and/or `tree_sample_prop` parameters in the `trees_cbh()` function. See (#10, @georgewoolsey).

- New: Adds the `leafr_for_ladderfuelsr()` which re-writes `leafR` package steps ([https://github.com/DRAAlmeida/leafR](https://github.com/DRAAlmeida/leafR)) to allow for an attribute input (e.g. "treeID") and removes the need to write individual tree point clouds to disk to extract LAD. The output of this function is used as input for the `ladderfuelsR` ([https://github.com/olgaviedma/LadderFuelsR](https://github.com/olgaviedma/LadderFuelsR)) steps as implemented in `ladderfuelsr_cbh()` in the present package.
- New: Adds the `polygon_attribute_to_las()` function to attach polygon attribute to point cloud
- Change: `ladderfuelsr_cbh()` now accepts input generated from `leafr_for_ladderfuelsr()` in addition to maintaining backwards capability
- Change: `trees_cbh()` incorporates the new process to extract CBH noted above and increases the default sample size
- Internal: R/check_las_data.R internal utility to check if input is a readable point cloud

# cloud2trees 0.5.0

Implements methods to estimate tree biomass (or crown biomass) using stand-level estimates (e.g. LANDFIRE raster estimates) distributed across the individual trees in the stand. See (#9, @georgewoolsey).

- New: Adds the function `trees_biomass()` to use an input tree list (e.g. as exported by `raster2trees()`) to estimate tree (or crown) biomass using one or many of the methods made available as listed in the documentation
- New: Adds the function `trees_biomass_cruz()` to use an input tree list (e.g. as exported by `raster2trees()`) to estimate tree crown biomass based on [Cruz et al. (2003)]((https://scholar.google.com/scholar?cluster=316241498622221569&oi=gsb&hl=en&as_sdt=0,5))
- New: Adds the function `trees_biomass_landfire()` to use an input tree list (e.g. as exported by `raster2trees()`) to estimate tree crown biomass based on [LANDFIRE's Forest Canopy Bulk Density (CBD)](https://landfire.gov/fuel/cbd) data
- New: Adds the function `trees_landfire_cbd()` to attach the raster cell value from [LANDFIRE's Forest Canopy Bulk Density (CBD)](https://landfire.gov/fuel/cbd) data to a tree list
- New: Adds the function `get_landfire()` to download the external [LANDFIRE's Forest Canopy Bulk Density (CBD)](https://landfire.gov/fuel/cbd) data
- Change: `get_data()` incorporates `get_landfire()`
- Change: `find_ext_data()` looks for LANDFIRE data
- Change: `cloud2trees()` incorporates the biomass process via `trees_biomass()` using the parameter `estimate_biomass_method`
- Internal: R/utils_biomass.R internal utility functions for distributing stand-level biomass to the tree-level
- Internal: R/check_spatial_points.R internal function to standardize process to check a tree list passed to the `trees_*()` functions

# cloud2trees 0.4.2

Updates point to raster matching functionality to allow for broader application and updates FIA Forest Type data to 30m resolution. See (#8, @georgewoolsey).

- Change: `get_foresttype()` now defaults to download 30m resolution data. If you previously ran `get_foresttype()` or `get_data()` and you want to use 30m data you will need to manually download the data again. See example below to update your local copy.
- New: Internal utility functions (see "utils_*.R") created to make point to raster data matching extendable to other data (e.g. LANDFIRE)
- Change: `trees_type()` implements new internal point to raster functions

To update your local copy of the FIA Forest Type data:

```r
cloud2trees::get_foresttype(force = T, res = 30)
```

# cloud2trees 0.4.1

- New: Makes the function `find_ext_data()` which attempts to locate the external data visible (was internal)
- Change: `cloud2raster()` writes the point cloud catalog file `raw_las_ctg_info.gpkg` to the `point_cloud_processing_delivery` directory

# cloud2trees 0.4.0

Implements a process to extract the height of the maximum crown diameter (HMD). See (#4, @georgewoolsey).

- New: Adds the function `trees_hmd()` to use input tree crown polygons (e.g. as exported by `raster2trees()`) to estimate tree HMD based on an input normalized point cloud
- New: Adds `simplify_multipolygon_crowns()` to simplify MULTIPOLYGON to POLYGON geometry in an `sf` class object
- Change: `cloud2trees()` incorporates the HMD process via `trees_hmd()`

# cloud2trees 0.3.1

Updates the processes that rely on the `lasR` package which implemented major revisions with breaking changes in [`lasR 0.13.0`](https://r-lidar.github.io/lasR/news/index.html#lasr-0130). The package `cloud2trees` now requires `lasR >= 0.13.1`. Even with these updates, users may continue to encounter some bugs in the near future. See (#3, @georgewoolsey).

To update execute:

```r
# update lasR
install.packages("lasR", repos = "https://r-lidar.r-universe.dev")
# update cloud2trees
remotes::install_github(repo = "georgewoolsey/cloud2trees", upgrade = F)
```

- Change: `lasr_*()` (internal functions) replace `LASlib` filters to conform with new `lasR` filtering (see [https://r-lidar.github.io/lasR/reference/filters.html](https://r-lidar.github.io/lasR/reference/filters.html))
- Change: `lasr_dtm_norm()` (internal function) replaces the `LASlib` filter `-keep_random_fraction` with `lasR::sampling_pixel()` to decimate the ground points (see [https://github.com/r-lidar/lasR/issues/102](https://github.com/r-lidar/lasR/issues/102))

# cloud2trees 0.3.0

Implements a process to extract the USDA Forest Inventory and Analysis (FIA) forest type group based on the [Forest Type Groups of the Continental United States](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d) data (Wilson 2023). See (#2, @georgewoolsey).

- New: Adds the function `trees_type()` to use the input tree list (e.g. as exported by `raster2trees()`) to attach species information using FIA codes
- New: Adds `get_data()` as all-in-one function that downloads all of the external data used by the package
- New: Adds `get_foresttype()` downloads the external data [Forest Type Groups of the Continental United States](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d)
- New: Adds `find_ext_data()` attempts to find the location of external data (as downloaded from `get_*()`)
- Change: `cloud2trees()` incorporates the CONUS forest type process via `trees_type()`

# cloud2trees 0.2.0

Integrates a process to extract the crown base height (CBH) using the workflow outlined in [Viedma et al. (2024)](https://doi.org/10.1111/2041-210X.14427) and using their package `LadderFuelsR` ([https://github.com/olgaviedma/LadderFuelsR](https://github.com/olgaviedma/LadderFuelsR)). See (#1, @georgewoolsey).

- New: Adds the function `trees_cbh()` to use input tree crown polygons (e.g. as exported by `raster2trees()`) to estimate tree CBH based on an input normalized point cloud
- New: Adds `ladderfuelsr_cbh()` as all-in-one function to extract CBH from a single tree, height normalized cloud using the functionality of the `LadderFuelsR` package
- Change: `cloud2trees()` incorporates the CBH process via `trees_cbh()`
- New: Adds error handling in `cloud2trees()` for stages after trees are extracted (e.g. during DBH extraction) the function will complete successfully and write the extracted trees data to the `point_cloud_processing_delivery` directory. An error message will be issued at the end (*ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! in: ...*) and state which portion did not complete successfully so that the user can attempt to re-run that portion with different settings (generally a `trees_*()` function)

# cloud2trees 0.1.0

- Open to public
