# internal functions to extract raster values at point locations

internal functions to extract raster values at point locations used by
[`trees_type()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_type.md)
and tree_biomass() for LANDFIRE rasters for the most part, these
functions take point and raster data as inputs the general process is
(see
[`trees_type()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_type.md)
for example):

- crop_raster_match_points()

- check for undesirable or NA values in the point_values from
  crop_raster_match_points()

- if undesirable values:

  - reclass_landfire_rast() or reclass_foresttype_rast()

  - agg_fill_rast_match_points(), see details above function definition

  - use the point_values from agg_fill_rast_match_points() if not null

  - otherwise, use the point_values from crop_raster_match_points()

Note for development: currently looking at entire bounding box of points
(or study boundary if defined) to fill NA raster cells this is the most
computationally intensive process, especially for large areas (\>100k
ha), fine resolution rasters (\<=30m) if the points are uniformly or
randomly dispersed across the entire AOI, the current process is likely
the best process if the points are clustered in groups across the AOI,
future development could attempt to group points into clusters and
iterate over the raster filling process for the cluster groups which may
allow for finer resolution raster data as input (instead of aggregating)
see examples

## Usage

``` r
crop_raster_match_points(
  points,
  rast,
  study_boundary = NA,
  max_search_dist_m = 1000
)
```

## Arguments

- points:

  sf.

- rast:

  SpatRaster.

- study_boundary:

  sf.

- max_search_dist_m:

  numeric.

## Examples

``` r
if (FALSE) { # \dontrun{
# !!!!!!!!!!!!!!!! this is a starter example for the author if decide to cluster points
# !!!!!!!!!!!!!!!! and apply raster match, fill, aggregate separately for clusters
# Sample data
data <- data.frame(
    x = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    y = c(2, 4, 3, 1, 5, 6, 7, 8, 9, 10)
)

# Calculate distances
distances <- dist(data)

# Perform hierarchical clustering
hc <- hclust(distances, method = "complete")

# Calculate WSS for different numbers of clusters
wss <- numeric(10)
for (k in 1:10) {
  wss[k] <- sum(kmeans(data, centers = k)$withinss)
}

# Plot the elbow
plot(1:10, wss, type = "b", xlab = "Number of Clusters (k)",
     ylab = "Within-Cluster Sum of Squares (WSS)")

# Find the optimal number of clusters (this is subjective)
# You might need to visually inspect the plot
optimal_k <- which.min(diff(wss, differences = 2))

# Cut the dendrogram based on the optimal number of clusters
clusters <- cutree(hc, k = optimal_k)
} # }
```
