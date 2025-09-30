# cloud2trees 0.8.0

# cloud2trees 0.7.2

Updates to use TreeMap 2022 ([https://doi.org/10.2737/RDS-2025-0032](https://doi.org/10.2737/RDS-2025-0032)). New data enables possible future development to include species classification modeling and CBH allometric prediction.

Users should update to this new data using: 

```r
get_treemap(force = T)
```

- Change: `trees_dbh()` uses whichever version of TreeMap is installed by detecting unique structure of data available. A warning is given if TreeMap 2022 is not downloaded. Also, improved error messages in function so that users can better diagnose issues.

# cloud2trees 0.7.1

Updates how `cloud2trees_to_lanl_trees()` writes data, includes unit tests for `cloud2trees_to_lanl_trees()` (by testing each processing step directly), adds description of file outputs and `cloud2trees_to_lanl_trees()` process to README, and more

- Change: `cloud2trees_to_lanl_trees()` was updated to format the inputs of the *fuellist* file and the *Cloud2Trees_TreeList.txt* file:
  + Rounded litter and grass bulk density to 3 digits in *fuellist* output file
  + Rounded litter and grass moisture to 2 digits in *fuellist* output file
  + Rounded litter and grass sizescale to 5 digits in *fuellist* output file
  + Rounded litter and grass depth to 2 digits in *fuellist* output file
  + Rounded all numeric data in the *Cloud2Trees_TreeList.txt* file to 4 digits
- Change: `cloud2trees_to_lanl_trees()` was updated to issue warnings/errors if the DTM extent does not fully encompass the study area extent, a requirement for many fire modeling tools
  + If there are any NA’s in the DTM and the `topofile` argument is set to "flat", then
    * Generate all files except the “topo.dat” file and issue a warning about the missing data
  + If there are any NA’s in the DTM and the `topofile` argument is set to "dtm", then
    * Halt execution and issue an error about the missing data
- Change: `cloud2trees()` and `raster2trees()` now check the `ws` argument (for ITD window function) prior to any data processing to avoid failure after some processing steps have already been completed
- Fix: `simplify_multipolygon_crowns()` now accounts for cases when a single tree has multiple polygon parts that are equal in area and also the largest part so that the function returns exactly the same number of records as the input data

# cloud2trees 0.7.0

Enable formatting of `cloud2trees()` data outputs for [LANL TREES](https://github.com/lanl/Trees/) program which does formatting for fire modeling

- New: Adds the function `cloud2trees_to_lanl_trees()` to use outputs from `cloud2trees()` to generate inputs for LANL TREES program
- New: all utility functions for working with tree list data given an AOI as well as preparing data for LANL TREES are in R/utils_aoi.R

# cloud2trees 0.6.9

Saving models used to estimate missing HMD and CBH values. Note, in the actual missing value estimation many RF models are estimated and model averaging is used. However, only the first estimated model is saved in this new export functionality which does not fully represent the process used to fill in missing values.

- Change: `trees_hmd()` now saves the model used to estimate missing HMD values to the path specified in `outfolder` if the parameter `estimate_missing_hmd = T`
- Change: `trees_cbh()` now saves the model used to estimate missing HMD values to the path specified in `outfolder` if the parameter `estimate_missing_hmd = T`
- Change: `cloud2trees()` automatically saves the model used to estimate missing HMD and CBH values to the `point_cloud_processing_delivery` directory

# cloud2trees 0.6.8

- Change: `cloud2raster()` now implements noise removal from the point cloud prior to the ground classification stage (in addition to after the stage) to minimize the influence of egregious outlier points on classification

# cloud2trees 0.6.7

- Fix: `itd_tuning()` now accounts for possible data transformation that occurs in `raster2trees()` on State Plane Coordinate System (SPCS) zone projections which use U.S. survey feet to express eastings and northings

# cloud2trees 0.6.6

Some point cloud data -- mostly USGS ALS acquisitions -- are made publicly available in State Plane Coordinate System (SPCS) zone projections which use U.S. survey feet to express eastings and northings (e.g. [EPSG:6430](https://epsg.io/6430)). Furthermore, some of these data utilize a coordinate system that combines two or more other coordinate systems, such as a horizontal and a vertical system by defining a Well-Known Text (WKT) using a Compound Coordinate System ("COMPD_CS"). For example, "NAD83(2011) / Colorado North (ftUS) + NAVD88 height - Geoid18 (ftUS)". None of the `cloud2trees` methods nor methods utilized by the program are designed to work with U.S. survey feet units. Furthermore, there is no simple way to transform data projected utilizing a compound coordinate system and transformations require manipulation of the full point cloud XYZ data ([see here](https://github.com/r-lidar/lidR/issues/372)). This update implements steps to look specifically for data projected using a horizontal projection in U.S. survey feet units, if detected, the program manipulates the full point cloud data to perform transformation. This transformation applies the [EPSG:5070](https://epsg.io/5070) projection, or "NAD83/Conus Albers", a coordinate reference system (CRS) with units in meters spanning the continental US (thus, from feet transformation will not work elsewhere). After transformation, the `cloud2trees` methods and methods utilized by the program work as designed on data projected using metric units. This workflow has been tested in two different SPCS zones with no issues identified and all unit tests are passing but issues may persist for novel point cloud data.

- New:  `raster2trees()`, and by extension `cloud2trees()`, implements steps to look specifically for data projected using a horizontal projection in U.S. survey feet units, if detected, the program manipulates the full point cloud data to perform transformation without any input from the user.
- New: all utility functions for manipulation of feet projections are in R/utils_projection.R

# cloud2trees 0.6.5

- Fix: `trees_biomass()` would refuse to overwrite values in "landfire_" and "cruz_" columns if those columns already existed in the data passed to the `tree_list`. This update now forces the overwrite of the data in these columns if already present in the data.

# cloud2trees 0.6.4

- Fix: `ladderfuelsr_cbh()` would fail when using the `las` parameter due to the lack of proper reference to the `pointsByZSlice()` function from the `leafR` package ([https://github.com/DRAAlmeida/leafR](https://github.com/DRAAlmeida/leafR)). This error is within the `leafR` package and as a workaround we define the global variable as `pointsByZSlice <<- leafR::pointsByZSlice`. Eventually, `ladderfuelsr_cbh()` needs to break the reliance on the `leafR` package.
- Fix: implements the internal function `as_character_safe()` defined in R/check_spatial_points.R to convert numeric columns to character which are meant to be used as a table identifier (i.e. as in "primary key") by ensuring that the number is not cast as character in scientific notation. To see the difference, compare `as.character(1000000000)` to `as_character_safe(1000000000)`.

# cloud2trees 0.6.3

- Fix: `trees_biomass()` (and the `trees_biomass_*()` functions) would return NA biomass values for trees that fell within a raster cell that *exactly* bordered the extent of the tree list. This change calculates the biomass for these border trees as if half of the raster cell (i.e. forest "stand") is within the study extent by updating the `calc_rast_cell_overlap()` utility function in R/utils_biomass.R

# cloud2trees 0.6.2

Several methods for attaching tree component metrics involve modelling missing values using a random forest model. The computational cost of random forests is driven by the repeated tree building process, which involves recursive partitioning, bootstrapping, and feature subset selection. When performed on XXL tree lists (e.g. 100k+) these operations result in a significant computational burden. Furthermore, large data sets can exceed the available RAM, leading to disk swapping, which significantly slows down the computation. To mitigate those memory problems when using `randomForest::randomForest()`, these updates implement data subsampling for large data sets to reduce memory consumption by randomly sampling a representative subset of the training data. This process is iterated with different subsets and the results are combined via model averaging. Model averaging is a technique for improving the robustness and accuracy of random forest models, especially when dealing with large data sets.

- Change: `trees_hmd()` now implements the random forest tuning and modelling techniques described above for large data sets
- Change: `trees_cbh()` now implements the random forest tuning and modelling techniques described above for large data sets
- Change: moves all utility functions for random forest models to R/utils_rf.R

# cloud2trees 0.6.1

- Change: same updates for `trees_hmd()` as in `0.6.0` which now allows for processing a list of files with crown polygon data such as the split crown polygon files "final_detected_crowns\*" written automatically for tree lists >250k by `raster2trees()` and `cloud2trees()`. This enables extracting HMD for XXL tree lists (e.g. 100k+) and effectively mitigates memory issues associated with these lists. The sampling is still done considering the full tree list.
- Change: moves all utility functions for `trees_hmd()` to R/utils_hmd.R

# cloud2trees 0.6.0

- Change: `trees_cbh()` now allows for processing a list of files with crown polygon data such as the split crown polygon files "final_detected_crowns\*" written automatically for tree lists >250k by `raster2trees()` and `cloud2trees()`. This enables extracting CBH for XXL tree lists (e.g. 100k+) and effectively mitigates memory issues associated with these lists. The sampling is still done considering the full tree list.
- Change: moves all utility functions for `trees_cbh()` to R/utils_cbh.R

The `trees_poly` parameter in `trees_cbh()` now accepts:

* `sf` class object with POLYGON geometry (see [sf::st_geometry_type()]). Recommended for smaller tree lists (e.g. <100k) that can fit in memory.
* character vector with the path to a single or multiple spatial files that can be read by [sf::st_read()] and have POLYGON geometry. Recommended for large tree lists (e.g. 100k+) that might cause memory issues.
* character with the path to a directory that has "final_detected_crowns\*" files from `cloud2trees()` or `raster2trees()`. Recommended for large tree lists (e.g. 100k+) that might cause memory issues.

# cloud2trees 0.5.9

- New: `itd_ws_functions()` makes a list of default functions that can be used for determining a variable window size for the detection of individual trees accessible.
- Change: `chunk_las_catalog()` would not allow for the processing of point clouds with "NA" CRS projection using `cloud2raster()`. Enabling the processing of point cloud data with "NA" projection.
- Fix: `trees_cbh()` had a potential memory leak that would cause memory full issues when processing XXL tree lists (e.g. 100k+) that would cause "Error: cannot allocate vector of size..." errors or just completely overwhelm the machine. Taking steps to remedy but might not ever be fully resolved.

# cloud2trees 0.5.8

Several methods for attaching tree component metrics involve modelling missing values using a random forest model. In `cloud2trees` these random forest models are tuned for each unique run using `randomForest::tuneRF()` which enables model tuning by searching for the optimal `mtry` parameter (the number of variables randomly sampled as candidates at each split) using a cross-validation approach. However, computational cost increases significantly with the number of observations as `randomForest::tuneRF()` performs cross-validation internally for each `mtry` value it tries. With large model training data (e.g. 100k+ observations), each of these cross-validation runs involves building and evaluating many random forest trees, making the process very time-consuming. This update adds the internal function `rf_tune_subsample()` to implement steps to mitigate very long run-times when tuning random forests models. `rf_tune_subsample()` mitigates very long run-times by: 

1) Reducing the `ntreeTry` parameter to a smaller value. Tuning will be less precise, but it will finish in a reasonable time. The `ntree` parameter is then increased for the final model. 
2) Subsampling uses a smaller, representative subsample of the data (e.g., 10-20% of the data) to find a good `mtry` value on the subsample.


- Change: `trees_hmd()` implements `rf_tune_subsample()` to mitigate very long run-times
- Change: `trees_cbh()` implements `rf_tune_subsample()` to mitigate very long run-times
- Change: `trees_dbh()` implements `rf_tune_subsample()` to mitigate very long run-times

# cloud2trees 0.5.7

- Fix: `trees_cbh()` and `trees_hmd()` would not be able to match based on `treeID` if the `treeID` column was a character that could also be cast as a numeric value (e.g. "1111111") due to the writing of results to disk storage rather than keeping everything in memory by `lidR::catalog_apply()` and then re-reading of results. This update re-casts the `treeID` in it's original data type.
- Fix: `trees_cbh()` and `trees_hmd()` would fail when the `trees_poly` data contained so many trees (e.g. >1M) due to exceeding the maximum allowed size by the `future` package default options as set via `lidR::catalog_apply()`. This update chunks up XXL tree lists into groups of 500k for processing to limit the size constraint errors.
- Change: `trees_hmd()` now includes the  `tree_sample_n` and `tree_sample_prop` parameters to give the option to limit the sample size. Prior to this change, HMD had to *attempt* to be extracted for all trees.
- Change: `cloud2trees()` now includes the  `hmd_tree_sample_n` and `hmd_tree_sample_prop` parameters to give the option to limit the HMD sample size. Prior to this change, HMD had to *attempt* to be extracted for all trees. The limit for HMD samples via `cloud2trees()` is set to 20,000.
- Change: `itd_tuning()` renames the default window size functions to more appropriately match the function shape. The defaults are: exponential (`exp_fn`; concave up) which was `nonlin_fn` in previous versions, linear (`lin_fn`) which did not change, and logarithmic (`log_fn`; concave down) which was `exp_fn` in previous versions.


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
