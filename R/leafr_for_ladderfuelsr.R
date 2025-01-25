#' @title re-writes `leafR` steps to allow for treeID as input for `ladderfuelsR`
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
#' @return Returns an data.frame
#'
#' @internal
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
tree_voxel_metrics = function(las, grain.size = 1, id = NULL){
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
    las <- lidR::filter_poi(las, treeID == id)
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
      , func = ~list(N = length(Z))
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
      treeID = dplyr::coalesce(id, as.numeric(NA))
      , n = dplyr::coalesce(n,0)
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
      total_pulses = sum(n,na.rm = T)
      , cumsum_pulses = cumsum(n)
      #Pulses out for each voxel
      , pulses_out = total_pulses - cumsum_pulses
      #The pulses.out of voxel 1 is the pulses.in of voxel 2 and so on...
      , pulses_in = dplyr::lag(pulses_out, n = 1, default=dplyr::first(pulses_out))
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
      , all_numeric = T
    )
    # map over treeID for tree_voxel_metrics
    voxel_metrics <- las@data %>%
      dplyr::filter(!is.na(treeID)) %>%
      dplyr::pull(treeID) %>%
      unique() %>%
      purrr::map(\(x) tree_voxel_metrics(
        las = las, grain.size = voxel_grain_size_m, id = x
      )) %>%
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
leafr_lad_profile = function(
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
    , col_names = c(cols2group, "n") %>% unique()
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
      n = sum(n, na.rm = T)
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
      , total_pulses = sum(n, na.rm=T)
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

