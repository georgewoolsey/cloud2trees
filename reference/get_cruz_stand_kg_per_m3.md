# internal functions to estimate tree biomass

internal functions to extract raster values at point locations used by
`trees_biomass*()` functions for the most part, these functions are used
to distribute raster cell (i.e. "stand") level estimates of fuel load to
individual trees

## Usage

``` r
get_cruz_stand_kg_per_m3(
  forest_type_group_code,
  basal_area_m2_per_ha,
  trees_per_ha
)
```

## Arguments

- forest_type_group_code:

  numeric. as extracted by
  [`trees_type()`](https://georgewoolsey.github.io/cloud2trees/reference/trees_type.md)

- basal_area_m2_per_ha:

  numeric.

- trees_per_ha:

  numeric.
