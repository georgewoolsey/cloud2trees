# Create CHM and apply pits and spikes filling via `lasR`

use `lasR` to create CHM and apply pits and spikes filling for raster
based on St-Onge 2008 (see reference).

## Usage

``` r
lasr_chm(
  chm_file_name,
  chm_res = 0.25,
  min_height_m = 2,
  max_height_m = 70,
  lap_sz = 3
)
```

## Arguments

- chm_file_name:

  string. Where to write the CHM.

- chm_res:

  numeric. The desired resolution of the CHM produced in meters.

- min_height_m:

  numeric. Set the minimum height (m) for individual tree detection

- max_height_m:

  numeric. Set the maximum height (m) for the canopy height model

- lap_sz:

  numeric. Size of the Laplacian filter kernel (integer value, in
  pixels) for
  [`lasR::pit_fill()`](https://rdrr.io/pkg/lasR/man/pit_fill.html)

## Value

A `lasR` pipeline

## References

<https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=81365288221f3ac34b51a82e2cfed8d58defb10e>
