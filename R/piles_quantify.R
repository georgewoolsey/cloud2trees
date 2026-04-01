#' @title Spectral filtering of candidate slash piles
#'
#' @description Calculates area and diameter from the input polygons
#' and height, and bulk volume using the CHM data within the input polygon footprints.
#' The function utilizes the input data as-is and does not height-filter the CHM nor
#' area-filter the polygons, for example.
#' Both input data sets should be projected with a CRS using meters as the horizontal unit of measurement,
#' while the CHM should have meters as the vertical units.
#'
#' @param sf_poly An sf object of polygons to be quantified.
#' @param chm_rast A SpatRaster object representing the Canopy Height Model.
#'
#' @return Returns a list of objects:
#'
#' * sf_data = spatial data frame with the input polygons and these columns added: `area_m2`, `diameter_m`, `volume_m3` , `max_height_m`
#' * area_rast = an area raster adjusted for variations in pixel area caused by the curvature of the Earth and the coordinate reference system
#' * volume_rast = a raster where values represent the cell area multiplied by the CHM height to create a volume raster
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#'  ###########################################
#'  # piles_quantify() test
#'  ###########################################
#'  # load a chm
#'  chm_fnm <- system.file(package = "cloud2trees", "extdata", "piles_chm.tif")
#'  my_chm <- terra::rast(chm_fnm)
#'  terra::plot(my_chm, axes = F)
#'  # load polygon data
#'  polys_fnm <- system.file(package = "cloud2trees", "extdata", "piles_poly.gpkg")
#'  my_polys <- sf::st_read(polys_fnm, quiet = T)
#'  my_polys
#'  my_polys %>%
#'    terra::vect() %>%
#'    terra::plot(border = "cyan", lwd = 2, col = NA, add = T)
#'  # piles_quantify() that
#'  piles_quantify_ans <- piles_quantify(
#'    sf_poly = my_polys
#'    , chm_rast = my_chm
#'  )
#'  # what?
#'  piles_quantify_ans
#'  # here are pile polygons
#'  dplyr::glimpse(piles_quantify_ans$sf_data)
#'  # define plot fn for each plot to have its own fill scale
#'  my_plot_fn <- function(data, fill_var, pal = "Blues", plot_title) {
#'    ggplot2::ggplot(data = data) +
#'      ggplot2::geom_sf(mapping=ggplot2::aes(fill = .data[[fill_var]])) +
#'      ggplot2::scale_fill_distiller(palette = pal, direction = 1) +
#'      ggplot2::labs(subtitle = plot_title, fill = fill_var) +
#'      ggplot2::theme_light() +
#'      ggplot2::theme(
#'        legend.position = "top"
#'        , axis.text = ggplot2::element_blank()
#'        , axis.ticks = ggplot2::element_blank()
#'        , grid.panel = ggplot2::element_blank()
#'        , plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "bold")
#'      )
#'  }
#'  # plot for each metric returned
#'  p1 <- my_plot_fn(piles_quantify_ans$sf_data, fill_var = "area_m2", pal = "Blues", plot_title = "Area")
#'  p2 <- my_plot_fn(piles_quantify_ans$sf_data, fill_var = "max_height_m", pal = "Oranges", plot_title = "Height")
#'  p3 <- my_plot_fn(piles_quantify_ans$sf_data, fill_var = "diameter_m", pal = "Purples", plot_title = "Diameter")
#'  p4 <- my_plot_fn(piles_quantify_ans$sf_data, fill_var = "volume_m3", pal = "Greys", plot_title = "Volume")
#'  # combine with patchwork
#'  library(patchwork)
#'  patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2, nrow = 2)
#' }
#'
piles_quantify <- function(
    sf_poly
    , chm_rast
    # , sf_id = NA
) {
  # check polygons
  if(!inherits(sf_poly, "sf")){stop("must include `sf` data object in 'sf_poly'")}
  if( !all(sf::st_is(sf_poly, type = c("POLYGON", "MULTIPOLYGON"))) ){
    stop(paste0(
      "`sf_poly` data must be an `sf` class object with POLYGON geometry (see [sf::st_geometry_type()])"
    ))
  }
  sf_poly <- sf_poly %>% dplyr::ungroup()

  # check raster
  # convert to SpatRaster if input is from 'raster' package
  if(
    inherits(chm_rast, "RasterStack")
    || inherits(chm_rast, "RasterBrick")
  ){
    chm_rast <- terra::rast(chm_rast)
  }else if(!inherits(chm_rast, "SpatRaster")){
    stop("Input 'chm_rast' must be a SpatRaster from the `terra` package")
  }
  chm_rast <- chm_rast %>% terra::subset(subset = 1)
  if(
    as.numeric(terra::global(chm_rast, fun = "isNA")) == terra::ncell(chm_rast)
    # || as.numeric(terra::global(chm_rast, fun = "isNA")) >= round(terra::ncell(chm_rast)*0.98)
  ){
    stop("Input 'chm_rast' has all missing values")
  }

  # # check id
  # if(!inherits(sf_id, "character")){
  #   # stop("must include 'sf_id' as the unique identifier")
  #   sf_poly <- sf_poly %>%
  #     dplyr::mutate(idxxxxx = dplyr::row_number())
  #   sf_id <- "idxxxxx"
  # }else{
  #   if( !any( stringr::str_equal(names(sf_poly), sf_id) ) ){
  #     stop(paste0("could not locate '",sf_id,"' in sf_poly"))
  #   }
  # }

  # check overlap
  # Returns TRUE if any part of the vector geometry intersects the raster extent
  if(
    !any(terra::is.related(
      x = sf_poly %>%
        sf::st_transform(terra::crs(chm_rast)) %>%
        terra::vect()
      , y = terra::ext(chm_rast)
      , relation = "intersects"
    ))
  ){
    stop("Input 'sf_poly' does not overlap with 'chm_rast'")
  }
  #################################
  # area, volume of each cell
  #################################
  area_rast_temp <- terra::cellSize(chm_rast)
  names(area_rast_temp) <- "area_m2"
  # area_rast_temp %>% terra::plot()
  # then, multiply area by the CHM (elevation) for each cell to get a raster with cell volumes
  vol_rast_temp <- area_rast_temp*chm_rast
  names(vol_rast_temp) <- "volume_m3"
  # vol_rast_temp %>% terra::plot()
  #################################
  # zonal stats
  #################################
  # sum area within each segment to get the total area
  area_df_temp <- terra::zonal(
      x = area_rast_temp
      , z = sf_poly %>%
        sf::st_transform(terra::crs(chm_rast)) %>%
        terra::vect()
      , fun = "sum", na.rm = T
    ) %>%
    setNames("area_m2") %>%
    dplyr::mutate(area_m2 = dplyr::na_if(area_m2, NaN))
  # area_df_temp %>% dplyr::glimpse()
  # sum volume within each segment to get the total volume
  vol_df_temp <- terra::zonal(
      x = vol_rast_temp
      , z = sf_poly %>%
        sf::st_transform(terra::crs(chm_rast)) %>%
        terra::vect()
      , fun = "sum", na.rm = T
    ) %>%
    setNames("volume_m3") %>%
    dplyr::mutate(volume_m3 = dplyr::na_if(volume_m3, NaN))
  # vol_df_temp %>% dplyr::glimpse()
  # max ht within each segment to get the max ht
  ht_df_temp <- terra::zonal(
      x = chm_rast
      , z = sf_poly %>%
        sf::st_transform(terra::crs(chm_rast)) %>%
        terra::vect()
      , fun = "max", na.rm = T
    ) %>%
    setNames("max_height_m") %>%
    dplyr::mutate(max_height_m = dplyr::na_if(max_height_m, NaN))
  #################################
  # attach to sf
  #################################
  if(
    !identical(
      nrow(sf_poly)
      , nrow(area_df_temp)
      , nrow(vol_df_temp)
      , nrow(ht_df_temp)
    )
  ){
    stop("unable to find data in raster for given polygons")
  }

  ret_dta <- sf_poly %>%
    dplyr::select( -dplyr::any_of(c(
      "hey_xxxxxxxxxx"
      , "area_m2"
      , "volume_m3"
      , "max_height_m"
      , "volume_per_area"
    ))) %>%
    dplyr::bind_cols(
      area_df_temp
      , vol_df_temp
      , ht_df_temp
    ) %>%
    dplyr::mutate(
      # use area of raster so we know where the volume came from
      volume_per_area = volume_m3/area_m2
      # now return area from the polygon
      , area_m2 = sf::st_area(.) %>% as.numeric()
      # adjust the volume to account for missing chm data
      , volume_m3 = area_m2*volume_per_area
    )
  # ret_dta <- sf_poly %>%
  #   purrr::reduce(
  #     list(sf_poly, area_df_temp, vol_df_temp, ht_df_temp)
  #     , dplyr::left_join
  #     , by = sf_id
  #   ) %>%
  #   dplyr::mutate(
  #     volume_per_area = volume_m3/area_m2
  #   )

  # calculate diameter
  # st_calculate_diameter in utils_projection.R
  ret_dta <- ret_dta %>%
    dplyr::select( -dplyr::any_of(c(
      "hey_xxxxxxxxxx"
      , "diameter_m"
    ))) %>%
    st_calculate_diameter()

  return(
    list(
      sf_data = ret_dta
      , area_rast = area_rast_temp
      , volume_rast = vol_rast_temp
    )
  )
}
