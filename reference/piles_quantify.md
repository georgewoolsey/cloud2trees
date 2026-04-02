# Spectral filtering of candidate slash piles

Calculates area and diameter from the input polygons and height, and
bulk volume using the CHM data within the input polygon footprints. The
function utilizes the input data as-is and does not height-filter the
CHM nor area-filter the polygons, for example. Both input data sets
should be projected with a CRS using meters as the horizontal unit of
measurement, while the CHM should have meters as the vertical units.

## Usage

``` r
piles_quantify(sf_poly, chm_rast)
```

## Arguments

- sf_poly:

  An sf object of polygons to be quantified.

- chm_rast:

  A SpatRaster object representing the Canopy Height Model.

## Value

Returns a list of objects:

- sf_data = spatial data frame with the input polygons and these columns
  added: `area_m2`, `diameter_m`, `volume_m3` , `max_height_m`,
  `pct_chm_coverarge`, `volume_per_area`

  - where, `pct_chm_coverarge` represents the proportion of the CHM
    cells within the polygon boundary that are not NA. If this value is
    very low (e.g. \<50%), volume and height estimates are potentially
    not reliable.

  - and, `volume_m3` is always adjusted using `volume_per_area` to
    represent the bulk volume based on the non-NA data within the
    polygon

- area_rast = an area raster adjusted for variations in pixel area
  caused by the curvature of the Earth and the coordinate reference
  system

- volume_rast = a raster where values represent the cell area multiplied
  by the CHM height to create a volume raster

## Examples

``` r
if (FALSE) { # \dontrun{
 ###########################################
 # piles_quantify() test
 ###########################################
 # load a chm
 chm_fnm <- system.file(package = "cloud2trees", "extdata", "piles_chm.tif")
 my_chm <- terra::rast(chm_fnm)
 terra::plot(my_chm, axes = F)
 # load polygon data
 polys_fnm <- system.file(package = "cloud2trees", "extdata", "piles_poly.gpkg")
 my_polys <- sf::st_read(polys_fnm, quiet = T)
 my_polys
 my_polys %>%
   terra::vect() %>%
   terra::plot(border = "cyan", lwd = 2, col = NA, add = T)
 # piles_quantify() that
 piles_quantify_ans <- piles_quantify(
   sf_poly = my_polys
   , chm_rast = my_chm
 )
 # what?
 piles_quantify_ans
 # here are pile polygons
 dplyr::glimpse(piles_quantify_ans$sf_data)
 # how was our CHM coverage within the piles?
 summary(piles_quantify_ans$sf_data$pct_chm_coverarge)
 # define plot fn for each added metric to have its own fill scale
 my_plot_fn <- function(data, fill_var, pal = "Blues", plot_title) {
   ggplot2::ggplot(data = data) +
     ggplot2::geom_sf(mapping=ggplot2::aes(fill = .data[[fill_var]])) +
     ggplot2::scale_fill_distiller(palette = pal, direction = 1) +
     ggplot2::labs(subtitle = plot_title, fill = fill_var) +
     ggplot2::theme_light() +
     ggplot2::theme(
       legend.position = "top"
       , axis.text = ggplot2::element_blank()
       , axis.ticks = ggplot2::element_blank()
       , plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "bold")
     )
 }
 # plot for each metric returned
 p1 <- my_plot_fn(
   piles_quantify_ans$sf_data
   , fill_var = "area_m2", pal = "Blues"
   , plot_title = "Area"
 )
 p2 <- my_plot_fn(
   piles_quantify_ans$sf_data
   , fill_var = "max_height_m", pal = "Oranges"
   , plot_title = "Height"
 )
 p3 <- my_plot_fn(
   piles_quantify_ans$sf_data
   , fill_var = "diameter_m", pal = "Purples"
   , plot_title = "Diameter"
 )
 p4 <- my_plot_fn(
   piles_quantify_ans$sf_data
   , fill_var = "volume_m3", pal = "Greys"
   , plot_title = "Volume"
 )
 # combine with patchwork
 library(patchwork)
 patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2, nrow = 2)
} # }
```
