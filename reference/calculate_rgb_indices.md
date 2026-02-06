# Calculate a suite of RGB-based spectral indices and color space conversions

This internal helper function transforms a standard three-band RGB
raster into a multi-layered stack of spectral indices and color
components. These layers provide the basis for the spectral threshold
voting system. While
[`piles_detect()`](https://georgewoolsey.github.io/cloud2trees/reference/piles_detect.md)
identifies candidates based on structure, this function provides the
spectral data used to calculate the "inrange_th_votes" during the
[`piles_spectral_filter()`](https://georgewoolsey.github.io/cloud2trees/reference/piles_spectral_filter.md)
process. By evaluating indices like ExGR and color components like
Lab_a, the framework can effectively distinguish non-photosynthetic
slash from green biomass.

## Usage

``` r
calculate_rgb_indices(rgb_rast, red_band_idx, green_band_idx, blue_band_idx)
```

## Arguments

- rgb_rast:

  A multi-band SpatRaster containing Red, Green, and Blue bands.

- red_band_idx:

  Integer. The index of the red band in the input raster.

- green_band_idx:

  Integer. The index of the green band in the input raster.

- blue_band_idx:

  Integer. The index of the blue band in the input raster.

## Value

A SpatRaster object containing the 13 calculated spectral and color
space layers.

## Details

The function generates the following layers:

- **grvi_layer**: Green Red Vegetation Index.

- **rgri_layer**: Red Green Ratio Index.

- **vdvi_layer**: Visible Band-Difference Vegetation Index.

- **rgbvi_layer**: Red Green Blue Vegetation Index.

- **exg_layer**: Excess Green Index.

- **exr_layer**: Excess Red Index.

- **exgr_layer**: Excess Green-minus-Red Index.

- **hsv_hue, hsv_saturation, hsv_brightness**: Components from the
  Hue-Saturation-Value color model.

- **Lab_L, Lab_a, Lab_b**: Components from the CIELAB color space,
  specifically utilized for their ability to isolate green-to-red
  spectral shifts.

## Examples

``` r
if (FALSE) { # \dontrun{
 # load rgb data
 rgb_fnm <- system.file(package = "cloud2trees", "extdata", "piles_rgb.tif")
 my_rgb <- terra::rast(rgb_fnm)
 my_rgb
 terra::plotRGB(my_rgb)
 # calculate_rgb_indices() that
 calculate_rgb_indices_ans <- calculate_rgb_indices(
   rgb_rast = my_rgb
   , red_band_idx = 1
   , green_band_idx = 2
   , blue_band_idx = 3
 )
 # what?
 terra::nlyr(calculate_rgb_indices_ans)
 names(calculate_rgb_indices_ans)
 # look at all these indices
 calculate_rgb_indices_ans %>%
   terra::plot(
     nc = 5
     , col = grDevices::gray.colors(111, start = 0, end = 1)
     , mar = c(0.2,0.2,1.6,0.2)
     , axes = FALSE
     , legend = F
   )
} # }
```
