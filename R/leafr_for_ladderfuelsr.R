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
#' @param group_treeID logical. should output be grouped by treeID? If `TRUE` (default), the attribute `treeID` must exist
#' in the `las` data and must be numeric.
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
  , group_treeID = TRUE
  , relative = FALSE
) {
  ###### re-write of leafR::lad.voxels
  lad_voxels <- leafr_lad_voxels(
      las
      , voxel_grain_size_m = voxel_grain_size_m
      , k = k
      , group_treeID = group_treeID
    )

  ###### re-write of leafR::lad.profile()
  if(nrow(lad_voxels)==0){return(NULL)}
  if(group_treeID==T){
    lad_profile <- leafr_lad_profile(
        lad_voxels = lad_voxels
        , attribute = "treeID"
        , relative = relative
      )
  }else{
    lad_profile <- leafr_lad_profile(
        lad_voxels = lad_voxels
        , attribute = NULL
        , relative = relative
      )
  }

  return(lad_profile)
}


#####################################################################
# intermediate function 1:
# filter las for one tree, aggregate with lidR::voxel_metrics
# attach treeID to return data frame
#####################################################################
tree_voxel_metrics <- function(
  las
  , grain.size = 1
  , id = NA # id should be numeric
){
  ## !!!!!! this is a re-write of leafR::lad.voxels()
  ## names roughly match what is in that function

  # force only one id
  id <- id[1]
  # filter the point cloud if there is a treeID
  if(
    (names(las@data) %>% stringr::str_equal("treeID") %>% any()) &&
    !is.na(id) &&
    !is.null(id)
  ){
    las <- lidR::filter_poi(las, as.character(treeID) == as.character(id))
  }

  # return nothing
  if (lidR::is.empty(las)) return(NULL)

  # force below ground to zero
  las@data$Z[las@data$Z < 0] <- 0
  # summary(las@data$Z)

  # take the floor of the Z values as done in leafR::pointsByZSlice
  las@data$Z <- floor(las@data$Z)

  # check the grain size
  # floor the grain size (resolution) and don't allow to go below 1m
  grain.size <- grain.size %>%
    as.numeric() %>%
    floor() %>%
    max(1, na.rm = T)

  # "vertical resolution (Delta Z or Dz) was fixed at 1 m" https://doi.org/10.3390/rs11010092
  dz <- 1

  ###################################################
  # get voxel data lidR::voxel_metrics()
  ###################################################
  ## when res = 1: z = 0 is 0-1 meters; z = 1 is 1-2 meters;...;z = 11 is 11-12 meters
    # lidR::voxel_metrics(
    #     las = las
    #     , func = ~list(N = length(Z))
    #     , res = grain.size
    #   ) %>%
    #   lidR::plot(
    #     color="N"
    #     , pal = harrypotter::hp_pal(option = "slytherin",begin=0.3)
    #     , size = grain.size, bg = "white", voxel = TRUE
    #   )

  ### can only do this for one tree at a time because the fn does not use the "attribute" parameter
  ### make this a function to do one tree at a time, return df with treeID attached, bind rows
  ### , do the rest of the steps with group_by(treeID)
  t.binneds <-
    lidR::voxel_metrics(
      las = las
      , func = ~list(pulses = length(Z))
      , res = c(grain.size,dz)
      # , attribute = "treeID" # this doesn't do anything and only works for crown_metrics()
    ) %>%
    dplyr::rename_with(tolower)
  # t.binneds$z %>% summary()
  # t.binneds %>% View()

  # "vertical resolution (Delta Z or Dz) was fixed at 1 m" https://doi.org/10.3390/rs11010092
  c_dz <- seq(from=0, to=max(t.binneds$z,na.rm = T), by = dz)

  # fill out the voxel data so that every x,y has the same z levels
  t.binneds <- dplyr::tibble(z=c_dz) %>%
    tidyr::crossing(t.binneds %>% dplyr::distinct(x,y)) %>%
    dplyr::left_join(t.binneds, by = dplyr::join_by(x,y,z)) %>%
    dplyr::mutate(
      treeID = dplyr::coalesce(as.character(id), as.character(NA))
      , pulses = dplyr::coalesce(pulses,0)
    ) %>%
    dplyr::relocate(treeID,x,y) %>%
    dplyr::arrange(treeID,x,y,z)

  # t.binneds$z %>% summary()
  # t.binneds %>% dplyr::mutate(dplyr::across(c(x,y),as.character)) %>% kableExtra::kbl() %>% kableExtra::kable_styling()

  return(t.binneds)
}

# # example
# tree_voxel_metrics(las, grain.size = 2, id = 41) %>% View()
# las@data %>%
#   dplyr::filter(!is.na(treeID)) %>%
#   dplyr::pull(treeID) %>%
#   unique() %>%
#   .[1:11] %>%
#   purrr::map(\(x) tree_voxel_metrics(
#     las = las, grain.size = 1, id = x
#   )) %>%
#   dplyr::bind_rows() %>%
#   dplyr::glimpse()

#####################################################################
# intermediate function 2:
# take response from tree_voxel_metrics
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
    cols2group <- c(attribute, "x","y")
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
    dplyr::mutate(dplyr::across(
      c("x","y","z")
      , as.numeric
    )) %>%
    ### this is essential...at a minimum, need to group by x,y
    dplyr::group_by(
      dplyr::across(
        dplyr::any_of(cols2group)
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
  , group_treeID = T
){
  if(!inherits(las, "LAS")){stop("must use an object of class `LAS`")}

  if(group_treeID){
    # check_df_cols_all_missing() in utils_biomass.r
    check_df_cols_all_missing(
      las@data
      , col_names = "treeID"
      , all_numeric = F
    )
    # map over treeID for tree_voxel_metrics
    voxel_metrics <- las@data %>%
      dplyr::filter(!is.na(treeID)) %>%
      dplyr::pull(treeID) %>%
      unique() %>%
      purrr::map(\(x)
        tree_voxel_metrics(
          las = las, grain.size = voxel_grain_size_m, id = x
        )
        , .progress = "extracting LAD (leaf area density)"
      ) %>%
      dplyr::bind_rows()
    # aggregate
    agg_df <- agg_lad_voxels(voxel_metrics, attribute="treeID") %>%
      dplyr::arrange(treeID,x,y,desc(z))
  }else{
    # tree_voxel_metrics
    voxel_metrics <- tree_voxel_metrics(las = las, grain.size = voxel_grain_size_m)
    # aggregate
    agg_df <- agg_lad_voxels(voxel_metrics) %>%
      dplyr::arrange(x,y,desc(z))
  }

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
# intermediate function 3:
# take response from leafr_lad_voxels
# and force it to native leafR::lad.voxels() return
#####################################################################
format_leafr_lad_voxels <- function(leafr_lad_voxels, id = NA) {
  if(
    (leafr_lad_voxels %>%
      dplyr::ungroup() %>%
      dplyr::summarise(n = sum(is.na(treeID)))) == nrow(leafr_lad_voxels)
  ){
    LAD_VOXELS = list(
      "LAD" = leafr_lad_voxels %>%
        dplyr::select(x,y,z,LAD) %>%
        dplyr::arrange(x,y,desc(z)) %>%
        tidyr::pivot_wider(
          names_from = z
          , values_from = LAD
          , names_prefix = "pulses_"
        ) %>%
        dplyr::select(-c(x,y)) %>%
        as.matrix()
      , "coordenates" = leafr_lad_voxels %>%
        dplyr::distinct(x,y) %>%
        dplyr::arrange(x,y) %>%
        as.matrix()
    )
  }else{

    LAD_VOXELS = list(
      "LAD" = leafr_lad_voxels %>%
        dplyr::filter(treeID==id) %>%
        dplyr::select(x,y,z,LAD) %>%
        dplyr::arrange(x,y,desc(z)) %>%
        tidyr::pivot_wider(
          names_from = z
          , values_from = LAD
          , names_prefix = "pulses_"
        ) %>%
        dplyr::select(-c(x,y)) %>%
        as.matrix()
      , "coordenates" = leafr_lad_voxels %>%
        dplyr::filter(treeID==id) %>%
        dplyr::distinct(x,y) %>%
        dplyr::arrange(x,y) %>%
        as.matrix()
    )
  }
  return(LAD_VOXELS)
}

# # # #
# leafr_lad_voxels(las %>% lidR::filter_poi(treeID==6),group_treeID = T) %>%
#   format_leafr_lad_voxels(id=6)
# # # lets compare...we won't be exact b/c we will get a different x,y grid
# f <- las %>%
#   lidR::filter_poi(treeID==6) %>%
#   lidR::writeLAS(file.path(tempdir(),"tree.las"))
#
# LAD_VOXELS_og <- leafR::lad.voxels(f)
# LAD_VOXELS_og

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
    cols2group <- c(attribute, "z")
    overall_cols2group <- c(attribute, "hey_xxxxxx") # hey_xxxxxx is a placeholder in case blank and works b/c using dplyr::any_of
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
    ungroup()
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
cols2group_to_id <- function(df, cols2group) {
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
  # return
  return(df)
}


#####################################################################
# intermediate function 50:
# aggregate points in voxels
# is similar to lidR::voxel_metrics()
#####################################################################
voxel_count_pulses <- function(
  las
  , horizontal_res = 1 # cannot go below 1
  , vertical_res = 1 # cannot go below 1
  , attribute = NULL # grouping attribute, e.g. "treeID"
  , force_z_gte0 = T # force Z<0 to zero
  , trim_voxels = "xyz"
  ## trim_voxels = "xyz" -> removes voxels where the xyz combination has zero pulses
  ### that is, returns only voxels that contain 1 or more points
  ## trim_voxels = "xy" -> removes voxels where the xy combination has zero pulses
  ## trim_voxels = "z" -> removes voxels where the z level has zero pulses
  ## trim_voxels = "none" -> returns all voxels in the bounding box of the attribute
) {
  # get the xyz bounding box for each tree
  voxel_df <- voxelize_las_to_bbox_df(
    las = las
    , horizontal_res = horizontal_res # cannot go below 1
    , vertical_res = vertical_res # cannot go below 1
    , attribute = attribute
    , full_z_range = T # should the return have the full z range even if no points?
    , force_z_gte0 = force_z_gte0 # force Z<0 to zero
  )
  if(is.null(voxel_df) || nrow(voxel_df)<1){return(NULL)}

  # force cols2group to something
  if(inherits(attribute,"character")){
    cols2group <- c(attribute) %>% unique() %>% stringr::str_squish()
  }else{
    ### make a dummy attribute
    ### this way we can continue to use the
    ### cols2group syntax below
    las@data$dummy_xxx <- 1
    voxel_df$dummy_xxx <- 1
    cols2group <- c("dummy_xxx")
  }

  ## in order to join by our cols2group we need to create a single column
  ## that uniquely identifies the columns listed in cols2group
  ## even if length(cols2group)>1
  ## cols2group_to_id() combines the values in cols2group
  ## and returns the data with a new column called "idxxx"
  voxel_df <- voxel_df %>%
    cols2group_to_id(cols2group = cols2group)
  las_data <-
    las@data %>%
    # remove rows where attribute is fully missing
    # !!! only do this if cols2group is not null
    dplyr::filter(
      !dplyr::if_all(
        .cols = dplyr::all_of(cols2group)
        , .fns = ~ dplyr::coalesce(as.numeric(as.factor(.x)), 0) == 0
      )
    ) %>%
    cols2group_to_id(cols2group = cols2group)

  # now we join the voxel data to the las points and count
  my_join <- dplyr::join_by(
    idxxx
    , x_from<=X, x_to>X # x_min <= x < x_max
    , y_from<=Y, y_to>Y # y_min <= y < y_max
    , z_from<=Z, z_to>Z # z_min <= z < z_max
  )
  # join and aggregate
  voxel_pulses_df <- voxel_df %>%
    dplyr::rename_with(
      .cols = c(x,y,z)
      , .fn = ~ paste0(.x, "_from", recycle0 = TRUE)
    ) %>%
    dplyr::mutate(
      x_to = x_from+horizontal_res
      , y_to = y_from+horizontal_res
      , z_to = z_from+vertical_res
    ) %>%
    ## potential memory issues for very large point clouds?
    dplyr::left_join(
      las_data %>%
        dplyr::select(dplyr::all_of(c(
          "idxxx"
          , "X","Y","Z"
        ))) %>%
        dplyr::mutate(pulses=1)
      , by = my_join # cols2group
    ) %>%
    dplyr::rename_with(
      .cols = c(x_from,y_from,z_from)
      , .fn = ~ stringr::str_remove_all(.x, "_from")
    ) %>%
    dplyr::group_by(dplyr::across(
      dplyr::all_of(c(
        cols2group
        , "x","y","z"
      ))
    )) %>%
    dplyr::summarise(
      pulses = sum(dplyr::coalesce(pulses,0))
    ) %>%
    dplyr::ungroup()

  # drop the highest xyz voxels if they don't have any pulses
  # these get introduced in voxelize_las_to_bbox_df to ensure full coverage
  voxel_pulses_df <- voxel_pulses_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(cols2group))) %>%
    dplyr::mutate(dplyr::across(
      c(x,y,z)
      , .fns = ~ max(ifelse(pulses>0,.x,as.numeric(NA)), na.rm = T)
      , .names = "{.col}_maxxx"
    )) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      x<=x_maxxx
      , y<=y_maxxx
      , z<=z_maxxx
    ) %>%
    dplyr::select(-tidyselect::ends_with("_maxxx"))

  # voxel_pulses_df %>%
  #   dplyr::relocate(x,y,z) %>%
  #   dplyr::mutate(treeID = treeID %>% as.factor() %>% as.numeric()) %>%
  #   dplyr::rename_with(toupper) %>%
  #   lidR::LAS(crs = lidR::st_crs(las)) %>%
  #   lidR::plot(
  #     color="PULSES", breaks = "kmeans"
  #     , size = 1, bg = "white", voxel = TRUE
  #   )

  ##### trim voxels
  ## trim_voxels = "xyz" -> removes voxels where the xyz combination has zero pulses
    ### that is, returns only voxels that contain 1 or more points
  ## trim_voxels = "xy" -> removes voxels where the xy combination has zero pulses
  ## trim_voxels = "z" -> removes voxels where the z level has zero pulses
  ## trim_voxels = "none" -> returns all voxels in the bounding box of the attribute
  trim_voxels <- dplyr::coalesce(trim_voxels,"") %>%
    .[1] %>%
    tolower() %>%
    stringr::str_replace_all("[^[:alnum:]]", "")
  if(trim_voxels == "xyz"){
    voxel_pulses_df <- voxel_pulses_df %>%
      dplyr::filter(dplyr::coalesce(pulses,0)>0)
  }else if(trim_voxels == "xy"){
    voxel_pulses_df <- voxel_pulses_df %>%
      dplyr::group_by(dplyr::across(
        dplyr::all_of(c(cols2group,"x","y"))
      ))
  }else if(trim_voxels == "z"){
    voxel_pulses_df <- voxel_pulses_df %>%
      dplyr::group_by(dplyr::across(
        dplyr::all_of(c(cols2group,"z"))
      ))
  } # otherwise no filter

  return(voxel_pulses_df)

  # voxel_pulses_df %>%
  #   dplyr::filter(dplyr::coalesce(pulses,0)>0) %>%
  #   dplyr::filter(treeID == voxel_pulses_df$treeID[1111]) %>%
  #   dplyr::relocate(x,y,z) %>%
  #   dplyr::mutate(treeID = treeID %>% as.factor() %>% as.numeric()) %>%
  #   dplyr::rename_with(toupper) %>%
  #   lidR::LAS(crs = lidR::st_crs(las)) %>%
  #   lidR::plot(
  #     color="PULSES", legend = T
  #     , size = 1, bg = "white", voxel = TRUE
  #   )
  # voxel_pulses_df %>%
  #   dplyr::filter(dplyr::coalesce(pulses,0)>0) %>%
  #   dplyr::filter(treeID == voxel_pulses_df$treeID[1]) %>%
  #   dplyr::arrange(x,y,z)
  #
  # lidR::voxel_metrics(
  #   las = las %>% lidR::filter_poi(treeID==voxel_pulses_df$treeID[1])
  #   , func = ~list(pulses = length(Z))
  #   , res = c(horizontal_res,vertical_res)
  # ) %>%
  # dplyr::rename_with(tolower) %>%
  # dplyr::arrange(x,y,z)
}
