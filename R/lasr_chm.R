#' @title Create CHM and apply pits and spikes filling via `lasR`
#'
#' @description
#' use `lasR` to create CHM and apply pits and spikes filling for raster based on St-Onge 2008 (see reference).
#'
#' @param chm_file_name string. Where to write the CHM.
#' @param chm_res numeric. The desired resolution of the CHM produced in meters.
#' @param min_height_m numeric. Set the minimum height (m) for individual tree detection
#' @param max_height_m numeric. Set the maximum height (m) for the canopy height model
#' @param lap_sz numeric. Size of the Laplacian filter kernel (integer value, in pixels) for [lasR::pit_fill()]
#'
#' @references
#' https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=81365288221f3ac34b51a82e2cfed8d58defb10e
#'
#' @return A `lasR` pipeline
#'
#' @keywords internal
#'
lasr_chm <- function(
 chm_file_name
 , chm_res = 0.25
 , min_height_m = 2
 , max_height_m = 70
 , lap_sz = 3
){
  # chm
    #set up chm pipeline step
    # operators = "max" is analogous to `lidR::rasterize_canopy(algorithm = p2r())`
    # for each pixel of the output raster the function attributes the height of the highest point found
    lasr_chm <- lasR::rasterize(
      res = chm_res
      , operators = "max"
      , filter = paste0(
        "-drop_class 2 9 18 -drop_z_below "
        , min_height_m
        , " -drop_z_above "
        , max_height_m
      )
      , ofile = ""
      # , ofile = paste0(config$chm_dir, "/*_chm.tif")
    )
  # Pits and spikes filling for raster with algorithm from St-Onge 2008 (see reference).
    # !!! testing with low point density lidar data revealed that the default pit_fill algorithm was too aggressive
    # ... decrease the lasR::pit_fill Size of the Laplacian filter kernel (integer value, in pixels) for low density point clouds
    lasr_chm_pitfill <- lasR::pit_fill(raster = lasr_chm, lap_size = lap_sz, ofile = chm_file_name)
  # pipeline
    pipeline <- lasr_chm + lasr_chm_pitfill
    return(pipeline)
}
