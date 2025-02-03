#' @title re-writes `leafR` steps to allow for treeID as input for `ladderfuelsR` or [ladderfuelsr_cbh()]
#'
#' @description
#' `leafr_for_ladderfuelsr()` is a re-write of [leafR::lad.voxels()] and [leafR::lad.profile()] to:
#'
#' * removes the requirement to use a file written to disk
#' * allows for the calculation of LAD by the `treeID` attribute so that don't have to pass individual tree point clouds
#' * updates to the use of the latest `lidR` functionality and removes the use of `sp` and `raster` functions
#' * updates the function to `tidy` data manipulation
#'
#' @param las an object of class LAS that has been height normalized.
#' @param voxel_grain_size_m numeric. horizontal resolution (suggested 1 meter for lad profiles and 10 meters for LAI maps). See `grain.size` in [leafR::lad.voxels()]
#' @param k numeric. coefficient to transform effective LAI to real LAI (k = 1; for effective LAI)
#' @param attribute character. The column name of the attribute to group the return LAD profile by. Default is "treeID". The attribute (whatever it is defined as)
#' must exist in the `las` data
#' @param min_pulses numeric. minimum number of pulses required to return a record by attribute. set to zero (default) to leave data unfiltered.
#' @param relative logical. produce lad profile by relative total LAI values. Indicate when using effective LAI value.
#' if set to TRUE, lad value will be relative_lad; otherwise, lad value will be mean_lad
#'
#' @references
#' * [https://doi.org/10.3390/rs11010092](https://doi.org/10.3390/rs11010092)
#' Almeida, D. R. A. D., Stark, S. C., Shao, G., Schietti, J., Nelson, B. W., Silva, C. A., ... & Brancalion, P. H. S. (2019). Optimizing the remote detection of tropical rainforest structure with airborne lidar: Leaf area profile sensitivity to pulse density and spatial sampling. Remote Sensing, 11(1), 92.
#' [https://github.com/DRAAlmeida/leafR](https://github.com/DRAAlmeida/leafR)
#'
#' @return Returns an data.frame which can have multiple treeIDs to use as input for [ladderfuelsr_cbh()]
#'
#' @examples
#'  \dontrun{
#'  # polygon data
#'  f <- system.file(package = "cloud2trees","extdata","crowns_poly.gpkg")
#'  trees_poly <- sf::st_read(f)
#'  # simplify polygons
#'  trees_poly <- simplify_multipolygon_crowns(trees_poly)
#'  # point cloud data
#'  lf <- system.file(package = "cloud2trees","extdata","norm_las","RMNP_017_2018_normalize.las")
#'  las <- lidR::readLAS(lf)
#'  las@data %>% dplyr::glimpse()
#'  # polygon_attribute_to_las to attach treeID to las
#'  las <- polygon_attribute_to_las(las, trees_poly, force_crs = T, attribute = "treeID")
#'  las@data %>% dplyr::glimpse()
#'  # get the lad profile for each treeID
#'  lad_profile <- leafr_for_ladderfuelsr(
#'      las
#'      , voxel_grain_size_m = 1
#'      , k = 1
#'      , group_treeID = T
#'      , relative = F
#'    )
#'  dplyr::glimpse(lad_profile)
#'  }
#' @export
#'
#'
leafr_for_ladderfuelsr <- function(
  las
  , voxel_grain_size_m = 1
  , k = 1
  , attribute = "treeID"
  , min_pulses = 0 # minimum number of pulses required to return record by attribute
  , relative = FALSE
) {
  ###### re-write of leafR::lad.voxels
  lad_voxels <- leafr_lad_voxels(
      las
      , voxel_grain_size_m = voxel_grain_size_m
      , k = k
      , attribute = attribute
      , min_pulses = min_pulses
    )
  if(is.null(lad_voxels) || nrow(lad_voxels)==0){return(NULL)}
  # lad_voxels %>% dplyr::glimpse()

  ###### re-write of leafR::lad.profile()
  lad_profile <- leafr_lad_profile(
    lad_voxels = lad_voxels
    , attribute = attribute
    , relative = relative
  )
  if(is.null(lad_profile) || nrow(lad_profile)==0){return(NULL)}
  # lad_profile %>% dplyr::glimpse()

  return(lad_profile)
}

#####################################################################
# intermediate function 2:
# take response from voxel_count_pulses
# aggregate and apply MacArthur-Horn equation to get LAD
#####################################################################
agg_lad_voxels = function(voxel_metrics, attribute = NULL, k = 1){
  # return nothing
  if(
    !inherits(voxel_metrics, "data.frame") ||
    (inherits(voxel_metrics, "data.frame") && nrow(voxel_metrics)<1)
  ){return(NULL)}
  # cols to group by
  if(inherits(attribute,"character")){
    cols2group <- c(attribute, "x","y") %>% unique() %>% stringr::str_squish()
  }else{
    cols2group <- c("x","y")
  }

  # check_df_cols_all_missing() in utils_biomass.r
  check_df_cols_all_missing(
    voxel_metrics
    , col_names = c(cols2group,"z") %>% unique()
    , all_numeric = F
  )

  # vertical resolution (Delta Z or Dz) was fixed at 1 m
  dz <- 1

  # aggregate and calculate LAD
  pulse_df <-
    voxel_metrics %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dplyr::across(
      c("x","y","z")
      , as.numeric
    )) %>%
    ### this is essential...at a minimum, need to group by x,y
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(cols2group)
      )
    ) %>%
    dplyr::arrange(x,y,desc(z)) %>% ### this is essential for the sky be first
    dplyr::mutate(
      total_pulses = sum(pulses,na.rm = T)
      , cumsum_pulses = cumsum(pulses)
      #Pulses out for each voxel
      , pulses_out = total_pulses - cumsum_pulses
      #The pulses.out of voxel 1 is the pulses.in of voxel 2 and so on...
      #In the highest voxel (closest to sky) the pulses in is the total pulses
      , pulses_in = dplyr::lag(pulses_out, n = 1, default=dplyr::first(total_pulses))
      # MacArthur-Horn eqquation
      # LAD = ln(S_bottom/S_top)*(1/(dz*K))
      #k value for LAD equation
      , LAD = (log(pulses_in/pulses_out) * 1/k * 1/dz) %>%
        ifelse(
          is.nan(.) | is.infinite(.)
          , as.numeric(NA)
          , .
        )
    ) %>%
    #remove the first 1 meter close to the ground (and the ground too)
    dplyr::filter(z>0) %>%
    dplyr::ungroup()

  # pulse_df %>% kableExtra::kbl() %>% kableExtra::kable_styling()
  return(pulse_df)
}

#####################################################################
# intermediate function 2.9:
# re-writes leafR::lad.voxels()
#####################################################################
leafr_lad_voxels <- function(
  las
  , voxel_grain_size_m = 1
  , k = 1
  , attribute = "treeID" # NULL
  , min_pulses = 0 # minimum number of pulses required to return record by attribute
){
  if(!inherits(las, "LAS")){stop("must use an object of class `LAS`")}
  # "vertical resolution (Delta Z or Dz) was fixed at 1 m" https://doi.org/10.3390/rs11010092
  dz <- 1
  # count the pulses
  voxel_metrics <- voxel_count_pulses(
    las
    , horizontal_res = voxel_grain_size_m # cannot go below 1
    , vertical_res = dz # cannot go below 1
    , attribute = attribute # grouping attribute, e.g. "treeID"
    , min_pulses = min_pulses # minimum number of pulses required to return record by attribute
    , force_z_gte0 = T # force Z<0 to zero
    , trim_voxels = "xy"
    ## trim_voxels = "xyz" -> removes voxels where the xyz combination has zero pulses
    ### that is, returns only voxels that contain 1 or more points
    ## trim_voxels = "xy" -> removes voxels where the xy combination has zero pulses
    ## trim_voxels = "z" -> removes voxels where the z level has zero pulses
    ## trim_voxels = "none" -> returns all voxels in the bounding box of the attribute
  )
  if(is.null(voxel_metrics) || nrow(voxel_metrics)<1){return(NULL)}
  # voxel_metrics %>% dplyr::glimpse()
  # aggregate
  agg_df <- agg_lad_voxels(voxel_metrics, attribute=attribute, k = k)
  if(is.null(agg_df) || nrow(agg_df)<1){return(NULL)}
  # agg_df %>% dplyr::glimpse()
  # return
  return(agg_df)

}

# leafr_lad_voxels(las,group_treeID = F) %>%
#   dplyr::filter(n>0) %>%
#   dplyr::rename_with(toupper) %>%
#   lidR::LAS() %>%
#   # purrr::pluck("data")
#   lidR::plot(
#         color="LAD"
#         , pal = harrypotter::hp_pal(option = "slytherin",begin=0.3)
#         , size = 1, bg = "white", voxel = TRUE
#       )

# leafr_lad_voxels(las,group_treeID = T) %>%
#   dplyr::filter(
#     n>0 &
#     treeID == (
#       las@data %>% dplyr::filter(!is.na(treeID)) %>%
#         dplyr::slice_head(n=1) %>% dplyr::pull(treeID)
#     )
#   ) %>%
#   dplyr::rename_with(toupper) %>%
#   lidR::LAS() %>%
#   # purrr::pluck("data")
#   lidR::plot(
#         color="LAD"
#         , pal = harrypotter::hp_pal(option = "slytherin",begin=0.3)
#         , size = 1, bg = "white", voxel = TRUE
#         , legend = T
#       )


#####################################################################
# intermediate function 4:
# re-writes leafR::lad.profile()
# to ingest output from leafr_lad_voxels() or agg_lad_voxels()
#####################################################################
leafr_lad_profile <- function(
  lad_voxels # return from agg_lad_voxels
  , attribute = NULL
  , relative = FALSE
){
  # return nothing if not df or no rows
  if(
    !inherits(lad_voxels, "data.frame") ||
    (inherits(lad_voxels, "data.frame") && nrow(lad_voxels)<1)
  ){return(NULL)}
  # cols to group by
  if(inherits(attribute,"character")){
    cols2group <- c(attribute, "z") %>% unique() %>% stringr::str_squish()
    overall_cols2group <- c(attribute, "hey_xxxxxx") %>%
       unique() %>% stringr::str_squish()
    # hey_xxxxxx is a placeholder in case blank and works b/c using dplyr::any_of
  }else{
    cols2group <- c("z")
    overall_cols2group <- c("hey_xxxxxx") # hey_xxxxxx is a placeholder in case blank and works b/c using dplyr::any_of
  }

  # check_df_cols_all_missing() in utils_biomass.r
  check_df_cols_all_missing(
    lad_voxels
    , col_names = c(cols2group, "pulses") %>% unique()
    , all_numeric = F
  )

  # average the height profiles (e.g. from Z= 2-3m) across the x,y to get one record per height profile
  lad_profile <- lad_voxels %>%
    ### this is essential...at a minimum, need to group by z
    dplyr::group_by(
      dplyr::across(
        dplyr::any_of(cols2group)
      )
    ) %>%
    dplyr::summarise(
      pulses = sum(pulses, na.rm = T)
      , mean_lad = mean(LAD, na.rm = T)
    ) %>%
    # calculate relative lad
    dplyr::group_by(
      dplyr::across(
        dplyr::any_of(overall_cols2group)
      )
    ) %>%
    dplyr::mutate(
      relative_lad = mean_lad/sum(mean_lad, na.rm=T)*100
      , total_pulses = sum(pulses, na.rm=T)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange((dplyr::across(dplyr::any_of(cols2group))))

  # # for some reason the original leafR::lad.profile() added 0.5 to the Z height
  #   # here are all the steps with the superlatives cut out
  #   # VOXELS_LAD$LAD is a matrix with rows of x,y and cols of lad_1_2m,...,lad_11_12m
  #   # VOXELS_LAD[[1]] is the same thing as VOXELS_LAD$LAD...so, yeah
  #   t.lad.profile = apply(VOXELS_LAD$LAD, 2, mean, na.rm = TRUE)
  #   max_height = ncol(VOXELS_LAD[[1]]) + .5
  #   t.lad.profile = data.frame(
  #     height = seq(1.5, max_height)
  #     , lad = t.lad.profile[length(t.lad.profile):1]
  #   )

  lad_profile <- lad_profile  %>%
    # # for some reason the original leafR::lad.profile() added 0.5 to the Z height
    dplyr::mutate(z = z+0.5) %>%
    dplyr::rename(height=z)
    # # for some reason the original leafR::lad.profile() keeps lad values that are NaN
    # # if we want to replace with NA, pipe the section below
    # dplyr::mutate(dplyr::across(
    #   .cols = tidyselect::ends_with("_lad")
    #   , .fns = ~ ifelse(
    #       is.nan(.x) | is.infinite(.x)
    #       , as.numeric(NA)
    #       , .x
    #     )
    # ))

  if(relative==T){
    lad_profile <- lad_profile  %>%
      dplyr::mutate(lad = relative_lad)
  }else{
    lad_profile <- lad_profile  %>%
      dplyr::mutate(lad = mean_lad)
  }

  # return
  return(lad_profile)

}
#####################################################################
# intermediate function 42:
## this is similar to lidR::voxelize_points() but allows for
## the grouping by attribute (e.g. treeID)
## and returns a data.frame with full bounding box coverage
## (can be converted to LAS) using something like:
## lidR::LAS(df, crs = st_crs(las))
#####################################################################
voxelize_las_to_bbox_df <- function(
  las
  , horizontal_res = 1 # cannot go below 1
  , vertical_res = 1 # cannot go below 1
  , attribute = NULL # grouping attribute, e.g. "treeID"
  , min_pulses = 0 # minimum number of pulses required to return record by attribute
  , full_z_range = T # should the return have the full z range even if no points?
  , force_z_gte0 = T # force Z<0 to zero
  ## full_z_range=T & force_z_gte0=T -> for every attribute z will go from [0, max(z)] by attribute
  ## full_z_range=F & force_z_gte0=T -> for every attribute z will go from [max(min(z),0), max(z)] by attribute
  ## full_z_range=F & force_z_gte0=F -> for every attribute z will go from [min(z), max(z)] by attribute
  ## full_z_range=T & force_z_gte0=F -> for every attribute z will go from [min(z@fulldata), max(z)] by attribute
) {
  if(!inherits(las, "LAS")){stop("must provide and object of class LAS")}
  # return nothing
  if (lidR::is.empty(las)) return(NULL)

  if(force_z_gte0==T){
    # force below ground to zero
    las@data$Z[las@data$Z < 0] <- 0
    # summary(las@data$Z)
  }

  # # take the floor of the Z values as done in leafR::pointsByZSlice
  # las@data$Z <- floor(las@data$Z)

  # check the grain size
  # floor the grain size (resolution) and don't allow to go below 1m
  horizontal_res <- horizontal_res %>%
    as.numeric() %>%
    floor() %>%
    max(1, na.rm = T)

  # check the grain size
  # floor the grain size (resolution) and don't allow to go below 1m
  vertical_res <- vertical_res %>%
    as.numeric() %>%
    floor() %>%
    max(1, na.rm = T)

  # # "vertical resolution (Delta Z or Dz) was fixed at 1 m" https://doi.org/10.3390/rs11010092
  # dz <- 1

  ## make a data.frame of the voxels by treeID
  # cols to group by
  if(
    inherits(attribute,"character")
  ){
    cols2group <- c(attribute) %>% unique() %>% stringr::str_squish()
    # # check_df_cols_all_missing() in utils_biomass.r
    # check_df_cols_all_missing(
    #   las@data
    #   , col_names = cols2group %>% unique()
    #   , all_numeric = F
    #   , check_vals_missing = T
    # )
  }else{
    ### make a dummy attribute
    ### this way we can continue to use the
    ### cols2group syntax below
    las@data$dummy_xxx <- 1
    cols2group <- c("dummy_xxx")
  }

  # prep the data and filter the data for this attribute
  las_data <-
    las@data %>%
    dplyr::rename_with(.fn = tolower, .cols = c(X,Y,Z)) %>%
    dplyr::select(dplyr::any_of(c(cols2group,"x","y","z"))) %>%
    # remove rows where attribute is fully missing
      # !!! only do this if cols2group is not null
    dplyr::filter(
      !dplyr::if_all(
        .cols = dplyr::all_of(cols2group)
        , .fns = ~ dplyr::coalesce(as.numeric(as.factor(.x)), 0) == 0
      )
    )

  # filter for min pulses
  if(dplyr::coalesce(as.numeric(min_pulses),0)>0){
    las_data <- las_data %>%
    dplyr::group_by(
      dplyr::across( dplyr::all_of(cols2group) )
    ) %>%
    dplyr::filter(dplyr::n()>as.numeric(min_pulses)) %>%
    dplyr::ungroup()
  }

  if(nrow(las_data)<1){return(NULL)}

  # create df unique by attribute with range of xyz values
  range_df <- las_data %>%
    dplyr::group_by(
      dplyr::across( dplyr::all_of(cols2group) )
    ) %>%
    dplyr::summarise(dplyr::across(
      .cols = c(x,y,z)
      , .fns = list(min = ~ min(.x,na.rm=T), max = ~ max(.x,na.rm=T))
    )) %>%
    dplyr::ungroup()
  if(nrow(range_df)<1){return(NULL)}
  # expand the grid by attribute because a row is unique by attribute
  grid_df <- range_df %>%
    dplyr::mutate(
      # x,y
      x_min = round(floor(x_min))
      , x_max = round(floor(x_max)+horizontal_res)
      , y_min = round(floor(y_min))
      , y_max = round(floor(y_max)+horizontal_res)
      # we are always starting at 0 ground for vertical resolution
      , z_min = dplyr::case_when(
        full_z_range==T & force_z_gte0==T ~ 0
        , full_z_range==T & force_z_gte0==F ~ dplyr::coalesce(min(las@data$Z, na.rm = T),0)
        , full_z_range==F & force_z_gte0==F ~ dplyr::coalesce(round(floor(z_min)),0)
        , full_z_range==F & force_z_gte0==T ~ dplyr::coalesce(round(floor(z_min)),0)
      )
      , z_max = round(floor(z_max)+vertical_res)
    ) %>%
    dplyr::rowwise() %>% # this is key
    dplyr::mutate(
      x = list(seq(from=x_min, to=x_max, by = horizontal_res))
      , y = list(seq(from=y_min, to=y_max, by = horizontal_res))
      , z = list(seq(from=z_min, to=z_max, by = vertical_res))
    ) %>%
    dplyr::select(-c(tidyselect::ends_with("_min"),tidyselect::ends_with("_max"))) %>%
    tidyr::unnest(cols = x) %>%
    tidyr::unnest(cols = y) %>%
    tidyr::unnest(cols = z) %>%
    dplyr::ungroup()
  if(nrow(grid_df)<1){return(NULL)}
  # remove dummy
  if(names(grid_df) %>% stringr::str_equal("dummy_xxx") %>% any()){
    grid_df <- grid_df %>% dplyr::select(-dummy_xxx)
  }
  return(grid_df)
}

# voxelize_las_to_bbox_df(las = las, attribute = c("treeID")) %>%
#   dplyr::glimpse()
# voxelize_las_to_bbox_df(las = las, attribute = "treeID") %>%
#   dplyr::relocate(x,y,z) %>%
#   dplyr::mutate(treeID = treeID %>% as.factor() %>% as.numeric()) %>%
#   dplyr::rename_with(toupper) %>%
#   lidR::LAS(crs = lidR::st_crs(las)) %>%
#   lidR::plot(
#     color="TREEID"
#     , size = 1, bg = "white", voxel = TRUE
#   )
# voxelize_las_to_bbox_df(las = las, attribute = "treeID", full_z_range = F) %>%
#   dplyr::relocate(x,y,z) %>%
#   dplyr::mutate(treeID = treeID %>% as.factor() %>% as.numeric()) %>%
#   dplyr::rename_with(toupper) %>%
#   lidR::LAS(crs = lidR::st_crs(las)) %>%
#   lidR::plot(
#     color="TREEID"
#     , size = 1, bg = "white", voxel = TRUE
#   )
# voxelize_las_to_bbox_df(las = las) %>%
#   dplyr::relocate(x,y,z) %>%
#   dplyr::rename_with(toupper) %>%
#   lidR::LAS(crs = lidR::st_crs(las)) %>%
#   lidR::plot(
#     color="Z"
#     , size = 1, bg = "white", voxel = TRUE
#   )
# voxelize_las_to_bbox_df(las = las, full_z_range = F) %>%
#   dplyr::relocate(x,y,z) %>%
#   dplyr::rename_with(toupper) %>%
#   lidR::LAS(crs = lidR::st_crs(las)) %>%
#   lidR::plot(
#     color="Z"
#     , size = 1, bg = "white", voxel = TRUE
#   )

#####################################################################
# intermediate function 48:
## in order to join by our cols2group we need to create a single column
## that uniquely identifies the columns listed in cols2group
## even if length(cols2group)>1
## if we were only joining our las to voxel data based on cols2group..
## ...we can use dplyr::left_join(voxel,las, by = cols2group)
## ...however, we are also joining on the range of x,y,z values
## ...and dplyr::join_by(cols2group, x_from<=X, x_to>X,..., z_from<=Z, z_to>Z)
## ...fails because the "cols2group" is a character list (non standard evaluation)
### we're going to do the same thing for both data so fn it
#####################################################################
cols2group_to_id <- function(df, cols2group, as_factor = F) {
  if(!inherits(df,"data.frame")){return(NULL)}
  if(!inherits(cols2group,"character")){return(df)}
  if(length(cols2group)<1){return(df)}
  ## 1. cast the cols2group as character
  df <-
    df %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dplyr::across(
      dplyr::all_of(cols2group)
      , .fns = ~ factor(.x) %>% as.character()
      , .names = "{.col}_fctxxx"
    ))
  ## 2. combine the cols2group columns into one id column
  df <- df %>%
    # create a combination of the cols2group as an ID
    dplyr::bind_cols(
      df %>%
      dplyr::select(tidyselect::ends_with("_fctxxx")) %>%
      tidyr::unite(
        col = "idxxx"
        , tidyselect::ends_with("_fctxxx")
        , sep = "_zzz_"
        , remove = T
        , na.rm = F
      )
    ) %>%
    dplyr::select(-tidyselect::ends_with("_fctxxx"))

  if(as_factor==T){
    df <- df %>% dplyr::mutate(idxxx = as.factor(idxxx))
  }
  # return
  return(df)
}

#####################################################################
# intermediate function 50.1:
# aggregate points in voxels
# is similar to lidR::voxel_metrics()
#####################################################################
voxel_count_pulses <- function(
  las
  , horizontal_res = 1 # cannot go below 1
  , vertical_res = 1 # cannot go below 1
  , attribute = NULL # grouping attribute, e.g. "treeID"
  , min_pulses = 0 # minimum number of pulses required to return record by attribute
  , force_z_gte0 = T # force Z<0 to zero
  , trim_voxels = "xyz"
  ## trim_voxels = "xyz" -> removes voxels where the xyz combination has zero pulses
    ### that is, returns only voxels that contain 1 or more points
  ## trim_voxels = "xy" -> removes voxels where the xy combination has zero pulses
    # this expands the data to all possible xyz combinations for which there is
    # a value in the xy even if there is no value in the z
    # so, a column from the ground up for all xy with points
  ## trim_voxels = "z" -> removes voxels where the z level has zero pulses
    # this expands the data to all possible xyz combinations for which there is
    # a value in the z even if there is no value in the xy
    # so, a row/bbox for all z levels with points
) {
  # las?
    if(!inherits(las, "LAS")){stop("must provide and object of class LAS")}
  # return nothing
    if (lidR::is.empty(las)) return(NULL)

    if(force_z_gte0==T){
      # force below ground to zero
      las@data$Z[las@data$Z < 0] <- 0
      # summary(las@data$Z)
    }

  # check the grain size
  # floor the grain size (resolution) and don't allow to go below 1m
    horizontal_res <- horizontal_res %>%
      as.numeric() %>%
      floor() %>%
      max(1, na.rm = T)

  # check the grain size
  # floor the grain size (resolution) and don't allow to go below 1m
    vertical_res <- vertical_res %>%
      as.numeric() %>%
      floor() %>%
      max(1, na.rm = T)

  # force cols2group to something
    if(inherits(attribute,"character")){
      cols2group <- c(attribute) %>% unique() %>% stringr::str_squish()
    }else{
      ### make a dummy attribute
      ### this way we can continue to use the
      ### cols2group syntax below
      las@data$dummy_xxx <- 1
      cols2group <- c("dummy_xxx")
    }

  # filter if missing all of the cols2group
    las_data <-
      las@data %>%
      # remove rows where attribute is fully missing
      # !!! only do this if cols2group is not null
      dplyr::filter(
        !dplyr::if_all(
          .cols = dplyr::all_of(cols2group)
          , .fns = ~ dplyr::coalesce(as.numeric(as.factor(.x)), 0) == 0
        )
      )

  # filter for min pulses
    if(dplyr::coalesce(as.numeric(min_pulses),0)>0){
      las_data <- las_data %>%
        dplyr::group_by(
          dplyr::across( dplyr::all_of(cols2group) )
        ) %>%
        dplyr::filter(dplyr::n()>as.numeric(min_pulses)) %>%
        dplyr::ungroup()
    }

    if(nrow(las_data)<1){return(NULL)}

  ## in order to join by our cols2group we need to create a single column
    ## that uniquely identifies the columns listed in cols2group
    ## even if length(cols2group)>1
    ## cols2group_to_id() combines the values in cols2group
    ## and returns the data with a new column called "idxxx"
    las_data <- las_data %>%
        cols2group_to_id(cols2group = cols2group, as_factor = T)
  # lookup table of grouping cols to id
    cols2group_to_id_df <- las_data %>%
      dplyr::select(dplyr::all_of(c(cols2group,"idxxx"))) %>%
      dplyr::distinct()

  ## aggregate the pulses by rounding of our xyz to the nearest multiple of resolution
    voxel_pulses_df <-
      las_data %>%
      dplyr::select(dplyr::all_of(c(
        "idxxx"
        , "X","Y","Z"
      ))) %>%
      dplyr::rename_with(tolower) %>%
      dplyr::group_by(idxxx) %>%
      # get the minimum or "starting" values
      dplyr::mutate(
        dplyr::across(
          c("x","y")
          , .fns = ~ min(.x,na.rm = T) %>%
            floor() %>%
            round_to_multiple(multiple = horizontal_res)
          , .names = "{.col}_min"
        )
        , z_min = min(z,na.rm = T) %>%
          floor() %>%
          round_to_multiple(multiple = vertical_res)
      ) %>%
      dplyr::ungroup() %>%
      # now group into voxels
      dplyr::mutate(
        x = round_to_multiple(x, multiple = horizontal_res, start = x_min)
        , y = round_to_multiple(y, multiple = horizontal_res, start = y_min)
        , z = round_to_multiple(floor(z), multiple = vertical_res, start = z_min) #leafR floors z
        , pulses=1
      ) %>%
      dplyr::group_by(idxxx,x,y,z) %>%
      dplyr::summarise(pulses = sum(pulses)) %>%
      dplyr::ungroup()

  # voxel_pulses_df %>% dplyr::glimpse()
  # voxel_pulses_df %>% summary()

  ### we are really expanding voxels here ;)
  ## trim_voxels = "xyz" -> removes voxels where the xyz combination has zero pulses
    ### that is, returns only voxels that contain 1 or more points
  ## trim_voxels = "xy" -> removes voxels where the xy combination has zero pulses
    # this expands the data to all possible xyz combinations for which there is
    # a value in the xy even if there is no value in the z
    # so, a column from the ground up for all xy with points
  ## trim_voxels = "z" -> removes voxels where the z level has zero pulses
    # this expands the data to all possible xyz combinations for which there is
    # a value in the z even if there is no value in the xy
    # so, a row/bbox for all z levels with points
  trim_voxels <- dplyr::coalesce(trim_voxels,"") %>%
    .[1] %>%
    tolower() %>%
    stringr::str_replace_all("[^[:alnum:]]", "")

  if(trim_voxels == "xyz"){
    # this is the same thing as no filter since started with las data
    # so a voxel has to have ate least one record
    voxel_pulses_df <- voxel_pulses_df %>%
      dplyr::filter(dplyr::coalesce(pulses,0)>0)
  }else if(trim_voxels == "xy"){
    voxel_pulses_df <-
      # get full range of z values by idxxx
      # from ground (or lowest point in the data) to highest point in group
      voxel_pulses_df %>%
        dplyr::ungroup() %>%
        dplyr::mutate(z_min = min(z,na.rm = T)) %>%
        dplyr::group_by(idxxx,z_min) %>%
        dplyr::summarise(z_max = max(z,na.rm = T)) %>%
        dplyr::ungroup() %>%
        dplyr::rowwise() %>% # this is key
        dplyr::mutate(
          z = list(seq(from=z_min, to=z_max, by = vertical_res))
        ) %>%
        dplyr::select(-c(tidyselect::ends_with("_min"),tidyselect::ends_with("_max"))) %>%
        tidyr::unnest(cols = z) %>%
        dplyr::ungroup() %>%
        # this expands the data to all possible xyz combinations for which there is
        # a value in the xy even if there is no value in the z
        # so, a column from the ground up for all xy with points
        dplyr::inner_join(
          # all x,y combinations by idxxx
          voxel_pulses_df %>%
            dplyr::distinct(idxxx,x,y)
          , by = "idxxx"
          , relationship = "many-to-many"
        ) %>%
        ## now join on voxel counts
        dplyr::left_join(
          voxel_pulses_df
          , by = dplyr::join_by(idxxx,x,y,z)
        ) %>%
        dplyr::mutate(pulses = dplyr::coalesce(pulses,0))
  }else if(trim_voxels == "z"){
    voxel_pulses_df <-
      # get full range of xy values (bbox) by idxxx
      voxel_pulses_df %>%
      dplyr::group_by(idxxx) %>%
      dplyr::summarise(dplyr::across(
        .cols = c(x,y)
        , .fns = list(min = ~ min(.x,na.rm=T), max = ~ max(.x,na.rm=T))
      )) %>%
      dplyr::ungroup() %>%
      dplyr::rowwise() %>% # this is key
      dplyr::mutate(
        x = list(seq(from=x_min, to=x_max, by = horizontal_res))
        , y = list(seq(from=y_min, to=y_max, by = horizontal_res))
      ) %>%
      dplyr::select(-c(tidyselect::ends_with("_min"),tidyselect::ends_with("_max"))) %>%
      tidyr::unnest(cols = x) %>%
      tidyr::unnest(cols = y) %>%
      dplyr::ungroup() %>%
      # this expands the data to all possible xyz combinations for which there is
      # a value in the z even if there is no value in the xy
      # so, a row/bbox for all z levels with points
      dplyr::inner_join(
        # all x,y combinations by idxxx
        voxel_pulses_df %>%
          dplyr::distinct(idxxx,z)
        , by = "idxxx"
        , relationship = "many-to-many"
      ) %>%
      ## now join on voxel counts
      dplyr::left_join(
        voxel_pulses_df
        , by = dplyr::join_by(idxxx,x,y,z)
      ) %>%
      dplyr::mutate(pulses = dplyr::coalesce(pulses,0))
  } # otherwise no filter which is same as all_voxels from lidR::voxel_metrics
  # attach id via lookup
  voxel_pulses_df <-
    voxel_pulses_df %>%
    dplyr::inner_join(
      cols2group_to_id_df
      , by = "idxxx"
      , relationship = "many-to-one"
    ) %>%
    dplyr::select(-c(idxxx)) %>%
    dplyr::relocate(dplyr::all_of(cols2group))
  # voxel_pulses_df %>% dplyr::glimpse()
  # voxel_pulses_df %>% summary()

  # voxel_pulses_df %>%
  #   dplyr::filter(pulses>0) %>%
  #   dplyr::filter(treeID == voxel_pulses_df$treeID[5555]) %>%
  #   dplyr::relocate(x,y,z) %>%
  #   dplyr::mutate(treeID = treeID %>% as.factor() %>% as.numeric()) %>%
  #   dplyr::rename_with(toupper) %>%
  #   lidR::LAS(crs = lidR::st_crs(las)) %>%
  #   lidR::plot(
  #     color="PULSES", legend = T
  #     , size = horizontal_res, bg = "white", voxel = TRUE
  #   )

  # return
  return(voxel_pulses_df)
}

#####################################################################
# intermediate function 53:
# round any number to the nearest multiple of a number
#####################################################################
round_to_multiple <- function(x, multiple, start=NULL) {
  remainder <- x %% multiple
  nearest <- ifelse(remainder < (multiple/2), x-remainder, x-remainder+multiple)
  y <- ifelse(nearest<dplyr::coalesce(start,nearest), start, nearest)
  return(y)
}
# round_to_multiple(c(0,1.5,5.3,11.0,24.1,25),2,4)
# round_to_multiple(c(0,1.5,5.3,11.0,24.1,25),2)
