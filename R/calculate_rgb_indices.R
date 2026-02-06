#' @title Calculate a suite of RGB-based spectral indices and color space conversions
#'
#' @description This internal helper function transforms a standard three-band
#' RGB raster into a multi-layered stack of spectral indices and color
#' components. These layers provide the basis for the spectral threshold
#' voting system. While [piles_detect()] identifies candidates based on
#' structure, this function provides the spectral data used to calculate
#' the "inrange_th_votes" during the [piles_spectral_filter()] process.
#' By evaluating indices like ExGR and color components like Lab_a,
#' the framework can effectively distinguish non-photosynthetic slash
#' from green biomass.
#'
#' @param rgb_rast A multi-band SpatRaster containing Red, Green, and Blue bands.
#' @param red_band_idx Integer. The index of the red band in the input raster.
#' @param green_band_idx Integer. The index of the green band in the input raster.
#' @param blue_band_idx Integer. The index of the blue band in the input raster.
#'
#' @details The function generates the following layers:
#' \itemize{
#'   \item \strong{grvi_layer}: Green Red Vegetation Index.
#'   \item \strong{rgri_layer}: Red Green Ratio Index.
#'   \item \strong{vdvi_layer}: Visible Band-Difference Vegetation Index.
#'   \item \strong{rgbvi_layer}: Red Green Blue Vegetation Index.
#'   \item \strong{exg_layer}: Excess Green Index.
#'   \item \strong{exr_layer}: Excess Red Index.
#'   \item \strong{exgr_layer}: Excess Green-minus-Red Index.
#'   \item \strong{hsv_hue, hsv_saturation, hsv_brightness}: Components from
#'   the Hue-Saturation-Value color model.
#'   \item \strong{Lab_L, Lab_a, Lab_b}: Components from the CIELAB color
#'   space, specifically utilized for their ability to isolate green-to-red
#'   spectral shifts.
#' }
#'
#' @return A SpatRaster object containing the 13 calculated spectral and color
#' space layers.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#'  # load rgb data
#'  rgb_fnm <- system.file(package = "cloud2trees", "extdata", "piles_rgb.tif")
#'  my_rgb <- terra::rast(rgb_fnm)
#'  my_rgb
#'  terra::plotRGB(my_rgb)
#'  # calculate_rgb_indices() that
#'  calculate_rgb_indices_ans <- calculate_rgb_indices(
#'    rgb_rast = my_rgb
#'    , red_band_idx = 1
#'    , green_band_idx = 2
#'    , blue_band_idx = 3
#'  )
#'  # what?
#'  terra::nlyr(calculate_rgb_indices_ans)
#'  names(calculate_rgb_indices_ans)
#'  # look at all these indices
#'  calculate_rgb_indices_ans %>%
#'    terra::plot(
#'      nc = 5
#'      , col = grDevices::gray.colors(111, start = 0, end = 1)
#'      , mar = c(0.2,0.2,1.6,0.2)
#'      , axes = FALSE
#'      , legend = F
#'    )
#' }
#'
calculate_rgb_indices <- function(rgb_rast, red_band_idx, green_band_idx, blue_band_idx) {
  # call individual index functions
  vdvi_layer <- spectral_index_vdvi(rgb_rast, red_band_idx, green_band_idx, blue_band_idx)
  grvi_layer <- spectral_index_grvi(rgb_rast, red_band_idx, green_band_idx)
  rgri_layer <- spectral_index_rgri(rgb_rast, red_band_idx, green_band_idx)
  rgbvi_layer <- spectral_index_rgbvi(rgb_rast, red_band_idx, green_band_idx, blue_band_idx)
  exg_layer <- spectral_index_exg(rgb_rast, red_band_idx, green_band_idx, blue_band_idx)
  exr_layer <- spectral_index_exr(rgb_rast, red_band_idx, green_band_idx, blue_band_idx)
  exgr_layer <- spectral_index_exgr(rgb_rast, red_band_idx, green_band_idx, blue_band_idx)
  # bi_layer <- spectral_index_bi(rgb_rast, red_band_idx, green_band_idx, blue_band_idx)
  # sat_layer <- spectral_index_saturation(rgb_rast, red_band_idx, green_band_idx, blue_band_idx)

  # color spaces: HSV
  hsv_rast <- spectral_rgb_to_hsv(rgb_rast, red_band_idx, green_band_idx, blue_band_idx)
  names(hsv_rast) <- paste0("hsv_", names(hsv_rast), recycle0 = T)
  hsv_hue <- hsv_rast$hsv_hue
  hsv_saturation <- hsv_rast$hsv_saturation
  hsv_brightness <- hsv_rast$hsv_brightness
  # color spaces: CEILab
  Lab_rast <- spectral_rgb_to_lab(rgb_rast, red_band_idx, green_band_idx, blue_band_idx)
  names(Lab_rast) <- paste0("Lab_", names(Lab_rast), recycle0 = T)
  Lab_L <- Lab_rast$Lab_L
  Lab_a <- Lab_rast$Lab_a
  Lab_b <- Lab_rast$Lab_b

  # stack all calculated indices into a single spatraster
  all_indices <- c(
    grvi_layer
    , rgri_layer
    , vdvi_layer
    , rgbvi_layer
    , exg_layer
    , exr_layer
    , exgr_layer
    # , bi_layer
    # , sat_layer
    , hsv_hue
    , hsv_saturation
    , hsv_brightness
    , Lab_L
    , Lab_a
    , Lab_b
  )

  return(all_indices)
}
###############################################################################
# check raster bands and extract rgb layers
###############################################################################
check_rgb_raster_bands <- function(rast, red_band_idx, green_band_idx, blue_band_idx) {
  # convert to SpatRaster if input is from 'raster' package
  if(
    inherits(rast, "RasterStack")
    || inherits(rast, "RasterBrick")
  ){
    rast <- terra::rast(rast)
  }else if(!inherits(rast, "SpatRaster")){
    stop("Input rgb raster must be a SpatRaster from the `terra` package")
  }

  # check if band indices are valid
  num_bands <- terra::nlyr(rast)
  # let 999999 be a cheat code
  cheat_code <- 999999

  if(
    ( red_band_idx!=cheat_code && red_band_idx > num_bands )
    || ( green_band_idx!=cheat_code && green_band_idx > num_bands )
    || ( blue_band_idx!=cheat_code &&  blue_band_idx > num_bands )
    || ( red_band_idx!=cheat_code && red_band_idx < 1 )
    || ( green_band_idx!=cheat_code && green_band_idx < 1 )
    || ( blue_band_idx!=cheat_code &&  blue_band_idx < 1 )
    || length(unique(c(red_band_idx,green_band_idx,blue_band_idx)))!=3
  ){
    stop("Invalid band index provided. Band indices must correspond to existing, unique layers in the rgb raster object.")
  }

  # extract bands
  if(red_band_idx!=cheat_code){ R <- rast[[red_band_idx]] }else{ R <- rast[[1]] }
  if(green_band_idx!=cheat_code){ G <- rast[[green_band_idx]] }else{ G <- rast[[1]] }
  if(blue_band_idx!=cheat_code){ B <- rast[[blue_band_idx]] }else{ B <- rast[[1]] }

  return(list(R = R, G = G, B = B))
}
# check_rgb_raster_bands(ortho_rast, red_band_idx=1, green_band_idx=3, blue_band_idx = 999999)
###############################################################################
# convert RGB to HSV color space
###############################################################################
spectral_rgb_to_hsv <- function(rast, red_band_idx, green_band_idx, blue_band_idx) {
  # check
  bands <- check_rgb_raster_bands(rast, red_band_idx, green_band_idx, blue_band_idx)
  R <- bands$R
  G <- bands$G
  B <- bands$B

  # convert to RGB
  rgb_rast <- terra::RGB(x = c(R,G,B), value = 1:3)

  # ?terra::colorize
  # convert to HSV (returns layers: h, s, v)
  hsv_rast <- terra::colorize(rgb_rast, to = "hsv")

  # names
  names(hsv_rast) <- c("hue", "saturation", "brightness")
  # scale only the hue layer to 360
  # hsv_rast$hue <- hsv_rast$hue*360

  return(hsv_rast)
}
###############################################################################
# convert RGB to CIELab color space
# ?grDevices::convertColor
###############################################################################
spectral_rgb_to_lab <- function(rast, red_band_idx, green_band_idx, blue_band_idx, max_dn = 255) {
  # check
  bands <- check_rgb_raster_bands(rast, red_band_idx, green_band_idx, blue_band_idx)
  R <- bands$R
  G <- bands$G
  B <- bands$B
  rgb_rast <- c(R,G,B)

  # values as a matrix (vectorized) normalized to max 255 val
  vals <- terra::values(rgb_rast) / max_dn

  # conversion on the matrix
  lab_vals <- grDevices::convertColor(vals, from = "sRGB", to = "Lab")

  # new SpatRaster using original str and replace values
  lab_rast <- terra::rast(rgb_rast, nlyrs = 3)
  terra::values(lab_rast) <- lab_vals

  names(lab_rast) <- c("L", "a", "b")

  return(lab_rast)
}
###############################################################################
###############################################################################
###############################################################################
# functions for individual RGB index
# functions for individual RGB index
# functions for individual RGB index
###############################################################################
###############################################################################
###############################################################################
# calculate green red vegetation index (GRVI)
# (G - R) / (G + R)
spectral_index_grvi <- function(rast, red_band_idx, green_band_idx) {
  bands <- check_rgb_raster_bands(rast, red_band_idx, green_band_idx, blue_band_idx = 999999) # blue_band_idx is dummy here
  R <- bands$R
  G <- bands$G

  grvi <- (G - R) / (G + R)
  names(grvi) <- "grvi"
  return(grvi)
}
# spectral_index_grvi(ortho_rast,red_band_idx = 1, green_band_idx = 2) %>% terra::plot()

# red green ratio index (RGRI)
# R/G
spectral_index_rgri <- function(rast, red_band_idx, green_band_idx) {
  bands <- check_rgb_raster_bands(rast, red_band_idx, green_band_idx, blue_band_idx = 999999) # blue_band_idx is dummy here
  R <- bands$R
  G <- bands$G

  rgri <- R/G
  names(rgri) <- "rgri"
  return(rgri)
}
# spectral_index_grvi(ortho_rast,red_band_idx = 1, green_band_idx = 2) %>% terra::plot()

# calculate visible band-difference vegetation index (VDVI)
# (2G - R - B) / (2G + R + B)
spectral_index_vdvi <- function(rast, red_band_idx, green_band_idx, blue_band_idx) {
  bands <- check_rgb_raster_bands(rast, red_band_idx, green_band_idx, blue_band_idx)
  R <- bands$R
  G <- bands$G
  B <- bands$B

  vdvi <- (2 * G - R - B) / (2 * G + R + B)
  names(vdvi) <- "vdvi"
  return(vdvi)
}
# spectral_index_vdvi(ortho_rast,red_band_idx = 1, green_band_idx = 2, blue_band_idx = 3) %>% terra::plot()

# calculate red green blue vegetation index (RGBVI)
# (G^2 - (B * R)) / (G^2 + (B * R))
spectral_index_rgbvi <- function(rast, red_band_idx, green_band_idx, blue_band_idx) {
  bands <- check_rgb_raster_bands(rast, red_band_idx, green_band_idx, blue_band_idx)
  R <- bands$R
  G <- bands$G
  B <- bands$B

  rgbvi <- (G^2 - (B * R)) / (G^2 + (B * R))
  names(rgbvi) <- "rgbvi"
  return(rgbvi)
}
# spectral_index_rgbvi(ortho_rast,red_band_idx = 1, green_band_idx = 2, blue_band_idx = 3) %>% terra::plot()

# calculate excess green (ExG)
# (2G - R - B) / (R + G + B) (using normalized RGB values)
# 2G - R - B (using raw values, then normalized by sum of R+G+B)
spectral_index_exg <- function(rast, red_band_idx, green_band_idx, blue_band_idx) {
  bands <- check_rgb_raster_bands(rast, red_band_idx, green_band_idx, blue_band_idx)
  R <- bands$R
  G <- bands$G
  B <- bands$B

  # Calculate normalized RGB values
  sum_rgb <- R + G + B
  r_norm <- R / sum_rgb
  g_norm <- G / sum_rgb
  b_norm <- B / sum_rgb

  exg <- (2 * g_norm - r_norm - b_norm)
  names(exg) <- "exg"
  return(exg)
}
# spectral_index_exg(ortho_rast,red_band_idx = 1, green_band_idx = 2, blue_band_idx = 3) %>% terra::plot()

# calculate brightness index (BI)
# sqrt((R^2 + G^2 + B^2) / 3)
spectral_index_bi <- function(rast, red_band_idx, green_band_idx, blue_band_idx) {
  bands <- check_rgb_raster_bands(rast, red_band_idx, green_band_idx, blue_band_idx)
  R <- bands$R
  G <- bands$G
  B <- bands$B

  bi <- sqrt((R^2 + G^2 + B^2) / 3)
  # normalize to 0-1 range
  max_brightness = max(
    terra::minmax(R)[2]
    , terra::minmax(G)[2]
    , terra::minmax(B)[2]
    , na.rm = T
  )

  bi <- bi/max_brightness
  names(bi) <- "bi"
  return(bi)
}
# spectral_index_bi(ortho_rast,red_band_idx = 1, green_band_idx = 2, blue_band_idx = 3) %>% terra::plot()

# calculate excessive red (ExR)
# 1.4r - g, where r and g are normalized RGB values.
spectral_index_exr <- function(rast, red_band_idx, green_band_idx, blue_band_idx){
  bands <- check_rgb_raster_bands(rast, red_band_idx, green_band_idx, blue_band_idx)
  R <- bands$R
  G <- bands$G
  B <- bands$B

  # Calculate normalized RGB values
  sum_rgb <- R + G + B
  r_norm <- R / sum_rgb
  g_norm <- G / sum_rgb

  exr <- (1.4 * r_norm - g_norm)
  names(exr) <- "exr"
  return(exr)
}

# calculate excess green-excess red (ExGR)
# 3g - 2.4r - b, where r, g, b are normalized RGB values.
# equivalent to ExG - ExR.
spectral_index_exgr <- function(rast, red_band_idx, green_band_idx, blue_band_idx) {
  bands <- check_rgb_raster_bands(rast, red_band_idx, green_band_idx, blue_band_idx)
  R <- bands$R
  G <- bands$G
  B <- bands$B

  # Calculate normalized RGB values
  sum_rgb <- R + G + B
  r_norm <- R / sum_rgb
  g_norm <- G / sum_rgb
  b_norm <- B / sum_rgb

  exgr <- (3 * g_norm - 2.4 * r_norm - b_norm)
  names(exgr) <- "exgr"
  return(exgr)
}

# calculate saturation (SAT)
# (max(R,G,B) - min(R,G,B)) / max(R,G,B)
spectral_index_saturation <- function(rast, red_band_idx, green_band_idx, blue_band_idx) {
  bands <- check_rgb_raster_bands(rast, red_band_idx, green_band_idx, blue_band_idx)
  R <- bands$R
  G <- bands$G
  B <- bands$B

  max_rgb <- max(R, G, B, na.rm = T)
  min_rgb <- min(R, G, B, na.rm = T)
  sat <- (max_rgb - min_rgb) / max_rgb
  names(sat) <- "sat"
  return(sat)
}
# spectral_index_saturation(ortho_rast,red_band_idx = 1, green_band_idx = 2, blue_band_idx = 3) %>% terra::plot()

