#' @title Estimate DBH for a tree list based on height
#'
#' @description
#' `trees_dbh()` uses the input tree list (e.g. as exported by [raster2trees()]) with the columns
#' `treeID`, `tree_x`, `tree_y`, and `tree_height_m` to estimate tree DBH.
#'
#' A regional model of height estimating DBH is determined by the process:
#'
#' * Use the TreeMap ([get_treemap()]) FIA plot data in the area of the tree list to estimate the height-DBH allometry relationship
#' * Use the height predicting DBH model built from the FIA data to predict DBH based on tree height in the tree list
#'
#' If training data is provided in `treels_dbh_locations` as returned from [treels_stem_dbh()]:
#'
#' * The regional model using FIA plot data is used to filter the DBH training data estimated from the point cloud
#' * The training data is used to estimate the height-DBH allometry relationship
#' * Use the height predicting DBH model built from the point cloud training data to predict DBH based on tree height in the tree list
#'
#' @param tree_list data.frame. A data frame with the columns `treeID`, `tree_x`, `tree_y`, and `tree_height_m`.
#' If an `sf` class object with POINT geometry (see [sf::st_geometry_type()]), the program will use the data "as-is" and only require the `treeID` and `tree_height_m` columns.
#' @param crs string. A crs string as returned from [sf::st_crs()] or the EPSG code of the x,y coordinates.
#' Defaults to the crs of the `tree_list` data if of class "sf".
#' @param study_boundary sf. The boundary of the study are to define the area of the regional model.
#' If no boundary given, regional model will be built from location of trees in the tree list.
#' @param dbh_model `r lifecycle::badge("deprecated")` Use the `dbh_model_regional` or `dbh_model_local` argument instead.
#' @param dbh_model_regional string. Set the model to use for regional dbh-height allometry based on FIA tree measurements.
#' Can be "cr" for the Chapman-Richards formula (default) or "power" for power function
#' @param dbh_model_local string. Set the model to use for local dbh-height allometry based on provided DBH training data in `treels_dbh_locations`.
#' Can be "rf" for random forest or "lin" for linear
#' @param treels_dbh_locations sf. Return from [treels_stem_dbh()].
#' Must also provide crown polygons (as returned from [raster2trees()]) in the `tree_list` data as an `sf` class object with POLYGON geometry (see [sf::st_geometry_type()])
#' If a valid file is provided, will make DBH predictions based on this training data instead of from the regional model from the FIA data
#' @param boundary_buffer numeric. Set the buffer (m) for the study area boundary to filter the FIA plot data based on TreeMap
#' @param input_treemap_dir directory where Treemap 2016 exists. Use [get_treemap()] first.
#' @param outfolder string. The path of a folder to write the model data to
#'
#' @references
#' * [https://doi.org/10.2737/RDS-2025-0032](https://doi.org/10.2737/RDS-2025-0032)
#' Houtman, Rachel M.; Leatherman, Lila S. T.; Zimmer, Scott N.; Housman, Ian W.; Shrestha, Abhinav; Shaw, John D.; Riley, Karin L. 2025. TreeMap 2022 CONUS: A tree-level model of the forests of the conterminous United States circa 2022. Fort Collins, CO: Forest Service Research Data Archive.
#'
#' * [https://doi.org/10.3390/f13122077](https://doi.org/10.3390/f13122077)
#' Tinkham et al. (2022). Modeling the missing DBHs: Influence of model form on UAV DBH characterization. Forests, 13(12), 2077.
#'
#' @return Returns a spatial data frame of individual trees.
#'
#' @examples
#'  \dontrun{
#'  library(tidyverse)
#'  # example tree list
#'  tl <- dplyr::tibble(
#'      treeID = c(1:21)
#'      , tree_x = rnorm(n=21, mean = 458064, sd = 11)
#'      , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
#'      , tree_height_m = exp(rgamma(n = 21, shape = (7/4)^2, rate = (4^2)/7))
#'    )
#'  # save our output somewhere (not required)
#'  outdir <- tempdir()
#'  # call the function
#'  tl_dbh <- trees_dbh(tree_list = tl, crs = "32613", outfolder = outdir)
#'  # what?
#'  tl_dbh %>% class()
#'  tl_dbh %>% dplyr::select(tidyselect::contains("dbh_cm")) %>% dplyr::glimpse()
#'  tl_dbh %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=dbh_cm))
#'  # what outputs did we get?
#'  list.files(outdir)
#'  # cloud2trees::trees_dbh() saved the FIA-measured trees used to train the allometric model
#'  read.csv(file.path(outdir, "regional_dbh_height_model_training_data.csv")) %>%
#'    summary()
#'  # cloud2trees::trees_dbh() saved the actual allometric model
#'  # let's load and review
#'  dbh_mod_temp <- readRDS(file.path(outdir, "regional_dbh_height_model.rds"))
#'  # what is this?
#'  dbh_mod_temp %>% class()
#'  # we can draw fit curves with probability bands using the tidybayes package
#'  library(tidybayes)
#'  # define our height range to predict over
#'  dplyr::tibble(tree_height_m = seq(from = 0, to = 25, by = 1)) %>%
#'    tidybayes::add_epred_draws(dbh_mod_temp, ndraws = 2000) %>%
#'    ggplot2::ggplot(ggplot2::aes(x = tree_height_m)) +
#'      tidybayes::stat_lineribbon(
#'        ggplot2::aes(y = .epred, color = "estimate")
#'        , .width = c(0.5,0.95)
#'        , lwd = 0.6
#'      ) +
#'      ggplot2::scale_fill_brewer(palette = "Oranges") +
#'      ggplot2::scale_color_manual(values = c("gray33")) +
#'      ggplot2::labs(x = "tree ht. (m)", y = "est. tree DBH (cm)", color = "") +
#'      ggplot2::scale_x_continuous(limits = c(0,NA), breaks = scales::extended_breaks(n=11)) +
#'      ggplot2::scale_y_continuous(limits = c(0,NA), breaks = scales::extended_breaks(n=11)) +
#'      ggplot2::theme_light()
#'  }
#' @export
#'
trees_dbh <- function(
  tree_list
  , crs = NA
  , study_boundary = NA
  , dbh_model_regional = "cr" # "power"
  , dbh_model_local = "lin"
  , treels_dbh_locations = NA
  , boundary_buffer = 50
  , input_treemap_dir = NULL
  , outfolder = tempdir()
){
  ####################################################################
  # check deprecated parameters
  ####################################################################
    calls <- names(sapply(match.call(), deparse))[-1]
    if(any("dbh_model" %in% calls)) {
        stop(
          "`dbh_model` deprecated. Use the `dbh_model_regional` or `dbh_model_local` argument instead."
        )
    }
  ####################################################################
  # check external data
  ####################################################################
    # find external data
    find_ext_data_ans <- find_ext_data(
      input_treemap_dir = input_treemap_dir
    )
    # if can't find external treemap data
    if(is.null(find_ext_data_ans$treemap_dir)){
      stop(paste0(
        "Treemap data has not been downloaded to package contents. Use `get_treemap()` first."
        , "\nIf you supplied a value to the `input_treemap_dir` parameter check that directory for data."
      ))
    }

    #treemap files
    treemap_data_finder_ans <- treemap_data_finder(find_ext_data_ans$treemap_dir)
    # treemap_data_finder_ans$which_treemap
    # treemap_data_finder_ans$treemap_trees
    # treemap_data_finder_ans$treemap_rast


  ##################################
  # ensure that tree height data exists
  ##################################
  f <- tree_list %>% names() %>% dplyr::coalesce("") # leaving in case anything below looks for this
  if(
    !(stringr::str_equal(f, "tree_height_m") %>% any())
  ){
    stop(paste0(
      "`tree_list` data must contain `tree_height_m` column to estimate DBH."
      , "\nRename the height column if it exists and ensure it is in meters."
    ))
  }
  tree_list <- tree_list %>%
    dplyr::mutate(tree_height_m = as.numeric(tree_height_m))

  ##################################
  # convert to spatial points data
  ##################################
  tree_tops <- check_spatial_points(tree_list, crs)
  if(sf::st_crs(tree_tops) %>% is.na()){
    stop(paste0(
      "Cannot make regional DBH-Height model with blank CRS."
      , "\n  ensure that the `tree_list` data has a CRS"
    ))
  }
  # get rid of columns we'll create
    tree_tops <- tree_tops %>%
      # throw in hey_xxxxxxxxxx to test it works if we include non-existant columns
      dplyr::select( -dplyr::any_of(c(
        "hey_xxxxxxxxxx"
        , "fia_est_dbh_cm"
        , "fia_est_dbh_cm_lower"
        , "fia_est_dbh_cm_upper"
        , "dbh_cm"
        , "is_training_data"
        , "dbh_m"
        , "radius_m"
        , "basal_area_m2"
        , "basal_area_ft2"
        , "ptcld_extracted_dbh_cm"
        , "ptcld_predicted_dbh_cm"
      )))

  ##################################
  # check for treels_dbh_locations
  ##################################
  if(dplyr::coalesce(nrow(treels_dbh_locations),0)>0){
    # check crown polygon data
    if( !(sf::st_is(tree_list, type = c("POLYGON", "MULTIPOLYGON")) %>% all()) ){
      stop(paste0(
        "`tree_list` data must be an `sf` class object with POLYGON geometry (see [sf::st_geometry_type()]) to estimate DBH from `treels_dbh_locations`"
        , "\nSee data returned from [raster2trees()] which is the intended data to pass to `tree_list`"
      ))
    }
    # check treels_dbh_locations
    if(!inherits(treels_dbh_locations, "sf")){
      stop(paste0(
        "`treels_dbh_locations` data must be an `sf` class object with POINT geometry (see [sf::st_geometry_type()]) to estimate DBH from `treels_dbh_locations`"
        , "\nSee data returned from [treels_stem_dbh()] which is the intended data to pass to `treels_dbh_locations`"
      ))
    }
    # check treels_dbh_locations
    if( !(sf::st_is(treels_dbh_locations, type = c("POINT", "MULTIPOINT")) %>% all()) ){
      stop(paste0(
        "`treels_dbh_locations` data must be an `sf` class object with POINT geometry (see [sf::st_geometry_type()]) to estimate DBH from `treels_dbh_locations`"
        , "\nSee data returned from [treels_stem_dbh()] which is the intended data to pass to `treels_dbh_locations`"
      ))
    }
    # check treels_dbh_locations
    if( !(names(treels_dbh_locations) %>% stringr::str_equal("dbh_cm") %>% any()) ){
      stop(paste0(
        "`treels_dbh_locations` data must have a column titled `dbh_cm` with numeric DBH values in cm."
        , "\nSee data returned from [treels_stem_dbh()] which is the intended data to pass to `treels_dbh_locations`"
      ))
    }
  }

  ##################################
  # define study boundary
  ##################################
  if(inherits(study_boundary, "sf") || inherits(study_boundary, "sfc")){
    buff <- study_boundary %>%
      sf::st_union() %>%
      sf::st_as_sf() %>%
      sf::st_transform(sf::st_crs(tree_tops)) %>%
      sf::st_buffer(boundary_buffer)
  }else{
    buff <- tree_tops %>%
      sf::st_bbox() %>%
      sf::st_as_sfc() %>%
      sf::st_buffer(boundary_buffer)
  }
  #check trees vs boundary
  nintersect <- tree_tops %>% sf::st_intersection(buff) %>% nrow() %>% dplyr::coalesce(0)
  if(nintersect==0){
    stop(paste0(
      "No trees in `tree_list` are within the `study_boundary`. DBH not estimated.\n"
      , " .... Check your data locations. If confident in tree locations, leave `study_boundary` as NA"
    ))
  }

  ####################################################################
  # DBH prediction using FIA data
  ####################################################################
    ###__________________________________________________###
    ### read in FIA data ###
    ###__________________________________________________###
    # can't make regional model if crs is blank
    proj_crs <- sf::st_crs(buff)$epsg
    proj_crs <- paste0("epsg:",tolower(proj_crs))
    if(proj_crs=="epsg:na"){
      stop(paste0(
        "Cannot make regional DBH-Height model with blank CRS."
        , "\nIf provided data in `study_boundary` ensure that it has a CRS"
        , "\nIf did not provide data in `study_boundary` ensure that the `tree_list` data has a CRS"
      ))
    }

    # read in treemap data
    # read in treemap (no memory is taken)
    treemap_rast <- terra::rast(treemap_data_finder_ans$treemap_rast)

    # check study boundary against the raster
    # the resulting matrix will have a true value where an intersection exists
    intersection_result <- terra::relate(
      x = buff %>%
        sf::st_union() %>%
        sf::st_transform(terra::crs(treemap_rast)) %>%
        terra::vect()
      , y = treemap_rast
      , relation = "intersects"
    )

    if(!any(intersection_result)) {
      stop(paste0(
        "The search area does not overlap with an area within CONUS. Cannot estimate DBH."
        , "\n.... If provided data in `study_boundary` ensure that it overlaps with CONUS"
        , "\n.... If did not provide data in `study_boundary` ensure that the `tree_list` data overlaps with CONUS"
      ))
    }

    ### filter treemap based on las...rast now in memory
    treemap_rast <- treemap_rast %>%
      terra::crop(
        buff %>%
          sf::st_union() %>%
          sf::st_transform(terra::crs(treemap_rast)) %>%
          terra::vect()
      ) %>%
      terra::subset(1)

    # ggplot(treemap_rast %>% as.data.frame(xy=T) %>% rename(f=3)) +
    #   geom_tile(aes(x=x,y=y,fill=as.factor(f))) +
    #   scale_fill_viridis_d(option = "turbo") +
    #   theme_light() + theme(legend.position = "none")

    # check for all NA's
    na_cells <- terra::global(treemap_rast, fun="isNA") %>% as.numeric()
    if(
      terra::ncell(treemap_rast)==dplyr::coalesce(na_cells,0)
    ){
      stop(paste0(
        "The search area does not overlap with a forested area within CONUS. Cannot estimate DBH."
        , "\n ... try expanding the `boundary_buffer` ?"
      ))
    }

    ### get weights for weighting each tree in the population models
    # treemap id = tm_id for linking to tabular data
    tm_id_weight_temp <- treemap_rast %>% ## works
      terra::values() %>%
      table() %>%
      dplyr::as_tibble() %>%
      dplyr::rename(tm_id=1,tree_weight=n) %>%
      dplyr::mutate(tm_id = as_character_safe(tm_id))
    # str(tm_id_weight_temp)

    ############################################################################
    ### get the TreeMap FIA tree list for only the plots included
    ############################################################################
    if(treemap_data_finder_ans$which_treemap==2022){
      treemap_cols <- c(
        "TM_ID"
        , "PLT_CN"
        , "SPECIES_SYMBOL"
        , "STATUSCD"
        , "DIA"
        , "HT"
        , "CR" # crown ratio
      )
    }else if(treemap_data_finder_ans$which_treemap==2016){
      treemap_cols <- c(
        "tm_id"
        , "CN"
        , "SPECIES_SYMBOL"
        , "STATUSCD"
        , "DIA"
        , "HT"
        , "CR" # crown ratio
      )
    }else{
      stop("unknown TreeMap vintage")
    }
    ### read it
    # treemap_data_finder_ans$treemap_trees
    treemap_trees_df <-
      readr::read_csv(
        treemap_data_finder_ans$treemap_trees
        , col_select = treemap_cols
        , progress = F
        , show_col_types = F
      ) %>%
      dplyr::rename_with(tolower)
    # rename cn to fit with original table str
    if(treemap_data_finder_ans$which_treemap==2022){
      treemap_trees_df <- treemap_trees_df %>% dplyr::rename(cn = plt_cn)
    }
    # clean it
    treemap_trees_df <-
      treemap_trees_df %>%
      dplyr::mutate(
        cn = as_character_safe(cn)
        , tm_id = as_character_safe(tm_id)
        , cr = ifelse(cr>100|cr<0,NA,cr)*0.01
      ) %>%
      dplyr::inner_join(
        tm_id_weight_temp
        , by = dplyr::join_by("tm_id")
      ) %>%
      dplyr::filter(
        # keep live trees only: 1=live;2=dead
        statuscd == 1
        & !is.na(dia)
        & !is.na(ht)
        & !is.na(tree_weight)
      ) %>%
      dplyr::mutate(
        dbh_cm = dia*2.54
        , tree_height_m = ht*0.3048
      ) %>%
      dplyr::select(-c(statuscd,dia,ht,cr)) # maybe we can use cr in a future version with and NSUR approach

    # check FIA model data
    if(dplyr::coalesce(nrow(treemap_trees_df),0)==0){
      stop(paste0(
        "Could not estimate DBH from Treemap data (see [`get_treemap()`]). Ensure that "
        , "\nIf provided data in `study_boundary` ensure that it is located in the continental US"
        , "\nIf did not provide data in `study_boundary` ensure that the `tree_list` data has is located in the continental US"
      ))
    }

    # save training data
    ### export tabular
      write.csv(
          treemap_trees_df %>% dplyr::mutate(which_treemap=treemap_data_finder_ans$which_treemap)
          , file.path(normalizePath(outfolder), "regional_dbh_height_model_training_data.csv")
          , row.names = F
        )

    ###__________________________________________________________###
    ### Regional model of DBH as predicted by height
    ### population model of dbh on height, non-linear
    ### used to filter sfm dbhs
    ###__________________________________________________________###
    if(stringr::str_squish( tolower(dbh_model_regional) )=="power"){
      # population model with no random effects (i.e. no group-level variation)
      # Define the non-linear model formula for DBH
      dbh_formula <- brms::bf(
        formula = dbh_cm|weights(tree_weight) ~ (b1 * tree_height_m) + tree_height_m^b2
        , b1 + b2 ~ 1
        , nl = TRUE # !! specify non-linear
      )

      # Define Priors
      dbh_priors <- c(
        brms::prior(normal(1, 2), nlpar = "b1")
        , brms::prior(normal(0, 2), nlpar = "b2")
      )

      # non-linear model form with Gamma distribution for strictly positive response variable dbh
      mod_nl_pop <- brms::brm(
        formula = dbh_formula
        , data = treemap_trees_df
        , prior = dbh_priors
        , family = brms::brmsfamily("Gamma")
        , iter = 6000, warmup = 3000, chains = 4
        , cores = lasR::half_cores()
        , file = paste0(normalizePath(outfolder), "/regional_dbh_height_model")
        , file_refit = "always"
      )
      # plot(mod_nl_pop)
      # summary(mod_nl_pop)
    }else{
      # dbh ~ asym * (1 - exp(-k * height))^p # Chapman-Richards non-linear formula
      # Define the non-linear Chapman-Richards formula for DBH
      dbh_formula <- brms::bf(
        dbh_cm|weights(tree_weight) ~ asym * (1 - exp(-k * tree_height_m))^p
        , asym ~ 1
        , k ~ 1
        , p ~ 1
        , nl = TRUE
      )

      # Define Priors
      # Note: Since we are using lognormal, the parameters asym, k, and p are
      # still interpreted on the original scale of the tree (e.g., inches or cm).
      dbh_priors <- c(
        # Asymptote: The max diameter. Adjust based on your units (e.g., inches vs cm).
        brms::prior(normal(60, 20), nlpar = "asym", lb = 0)
        # k: The rate of approach to the asymptote. Usually a small decimal.
        , brms::prior(normal(0.05, 0.02), nlpar = "k", lb = 0)
        # p: The shape parameter. p > 1 creates the initial acceleration.
        , brms::prior(normal(2, 0.5), nlpar = "p", lb = 0)
      )

      # non-linear model form with lognormal distribution for strictly positive response variable dbh
        # the lognormal family is appropriate for allometric data since
        # as trees get larger, the variance in their diameter typically increases.
        # the lognormal distribution naturally accounts for this because it assumes the
        # error is multiplicative rather than additive, and it ensures that predicted DBH values are always strictly positive.
      mod_nl_pop <- brms::brm(
        formula = dbh_formula
        , data = treemap_trees_df
        , prior = dbh_priors
        , family = brms::lognormal()
        , iter = 6000, warmup = 3000, chains = 4
        , cores = lasR::half_cores()
        , file = paste0(normalizePath(outfolder), "/regional_dbh_height_model")
        , file_refit = "always"
        , control = list(adapt_delta = 0.98)
      )
      # plot(mod_nl_pop)
      # summary(mod_nl_pop)
    }

    ## write out model estimates to tabular file
    #### extract posterior draws to a df
    brms::as_draws_df(
      mod_nl_pop
      , variable = c("^b_", "shape")
      , regex = TRUE
    ) %>%
      # quick way to get a table of summary statistics and diagnostics
      posterior::summarize_draws(
        "mean", "median", "sd"
        ,  ~quantile(.x, probs = c(
          0.05, 0.95
          , 0.025, 0.975
        ))
        , "rhat"
      ) %>%
      dplyr::mutate(
        variable = stringr::str_remove_all(variable, "_Intercept")
        , formula = summary(mod_nl_pop)$formula %>%
          as.character() %>%
          .[1]
      ) %>%
      write.csv(
        paste0(normalizePath(outfolder), "/regional_dbh_height_model_estimates.csv")
        , row.names = F
      )


    ### obtain model predictions over range
    # range of x var to predict
    height_range <- dplyr::tibble(
      tree_height_m = seq(
        from = 0
        , to = 120 # tallest tree in the world
        , by = 0.1 # by 0.1 m increments
      )
    )
    # predict and put estimates in a data frame
    pred_mod_nl_pop_temp <- predict(
      mod_nl_pop
      , newdata = height_range
      , probs = c(.05, .95)
    ) %>%
      dplyr::as_tibble() %>%
      dplyr::rename(
        lower_b = 3, upper_b = 4
      ) %>%
      dplyr::rename_with(tolower) %>%
      dplyr::select(-c(est.error)) %>%
      dplyr::bind_cols(height_range) %>%
      dplyr::rename(
        tree_height_m_tnth=tree_height_m
        , fia_est_dbh_cm = estimate
        , fia_est_dbh_cm_lower = lower_b
        , fia_est_dbh_cm_upper = upper_b
      ) %>%
      dplyr::mutate(tree_height_m_tnth=as_character_safe(tree_height_m_tnth)) %>%
      dplyr::relocate(tree_height_m_tnth)
    # str(pred_mod_nl_pop_temp)

    # save predictions for reading later
    write.csv(
      pred_mod_nl_pop_temp
      , file = paste0(normalizePath(outfolder), "/regional_dbh_height_model_predictions.csv")
      , row.names = F
    )

    # attach to treelist
    tree_tops <- tree_tops %>%
      # join with model predictions at 0.1 m height intervals
        dplyr::mutate(
          tree_height_m_tnth = round(as.numeric(tree_height_m),1) %>% as_character_safe()
        ) %>%
        dplyr::left_join(
          pred_mod_nl_pop_temp
          , by = dplyr::join_by(tree_height_m_tnth)
        ) %>%
        dplyr::select(-tree_height_m_tnth) %>%
        dplyr::mutate(dbh_cm = fia_est_dbh_cm)

  ####################################################################
  # Model DBH using `treels_dbh_locations`
  ####################################################################
    if(dplyr::coalesce(nrow(treels_dbh_locations),0)>0){
      ###________________________________________________________###
      ### Join the Top down crowns with the stem location points ###
      ###________________________________________________________###
        ### Join the top down crowns with the stem location points
        ## !! Note that one crown can have multiple stems within its bounds
      crowns_sf_joined_stems_temp <- tree_list %>%
        sf::st_join(
          treels_dbh_locations %>%
            # rename all columns to have "stem" prefix
            dplyr::rename_with(
              .fn = ~ paste0("stem_",.x,recycle0 = T)
              , .cols = tidyselect::everything()[
                -dplyr::any_of(
                  c(tidyselect::starts_with("stem_"),"stem_x", "stem_y","geom","geometry")
                )
              ]
            )
        )
      if(
        max(grepl("tree_x", names(crowns_sf_joined_stems_temp)))==0
        | max(grepl("tree_y", names(crowns_sf_joined_stems_temp)))==0
      ){ # doesn't contain x,y
        tree_tops <- crowns_sf_joined_stems_temp %>%
          sf::st_centroid() %>%
          dplyr::mutate(
            tree_x = sf::st_coordinates(.)[,1]
            , tree_y = sf::st_coordinates(.)[,2]
          )
      }
      ###________________________________________________________###
      ## Filter the SfM DBHs using the FIA model
      ###________________________________________________________###
        # 1) Predict an expected DBH value for each tree based on the FIA model
        # 2) Remove stem DBH estimates that are outside the 90% prediction bounds
        # 3) Select the stem DBH estimate that is closest to the predicted DBH value (from 1) if multiple stems are within the bounds of one crown
        # 4) Use the SfM-detected stems remaining after this filtering workflow as the training data in the local DBH to height allometric relationship model
        # attach allometric data to CHM derived trees and canopy data
      crowns_sf_joined_stems_temp <- crowns_sf_joined_stems_temp %>%
        # join with model predictions at 0.1 m height intervals
        dplyr::mutate(
          tree_height_m_tnth = round(as.numeric(tree_height_m),1) %>% as_character_safe()
        ) %>%
        dplyr::inner_join(
          pred_mod_nl_pop_temp
          , by = dplyr::join_by(tree_height_m_tnth)
        ) %>%
        dplyr::select(-tree_height_m_tnth) %>%
        dplyr::mutate(
          stem_dbh_cm = as.numeric(stem_dbh_cm)
          , fia_est_dbh_pct_diff = abs(stem_dbh_cm-fia_est_dbh_cm)/fia_est_dbh_cm
        )
        # what is the estimated difference
        # summary(crowns_sf_joined_stems_temp$fia_est_dbh_pct_diff)
        # crowns_sf_joined_stems_temp %>% dplyr::glimpse()
        # crowns_sf_joined_stems_temp %>% dplyr::filter(!is.na(stem_dbh_cm)) %>% nrow()

      ### build training data set by filtering stems
        dbh_training_data_temp <- crowns_sf_joined_stems_temp %>%
          sf::st_drop_geometry() %>%
          dplyr::filter(
            !is.na(stem_dbh_cm)
            & stem_dbh_cm > 0
            & stem_dbh_cm >= fia_est_dbh_cm_lower
            & stem_dbh_cm <= fia_est_dbh_cm_upper
          )
        if(nrow(dbh_training_data_temp)>0){
          # filter it again
          dbh_training_data_temp <- dbh_training_data_temp %>%
            dplyr::group_by(treeID) %>%
            # select the minimum difference to regional dbh estimate
            dplyr::filter(
              fia_est_dbh_pct_diff==min(fia_est_dbh_pct_diff, na.rm = T)
            ) %>%
            # just take one if same dbh
            dplyr::filter(
              dplyr::row_number()==1
            ) %>%
            dplyr::ungroup() %>%
            dplyr::select(c(
              treeID # id
              , stem_dbh_cm # y
              # x vars
              , tree_height_m
              , tree_x
              , tree_y
              # , crown_area_m2
              # , tidyselect::starts_with("min_crown_ht_m")
              # , tidyselect::starts_with("comp_")
            ))
        }
      ###__________________________________________________________###
      ### Build regional model to estimate missing DBHs using SfM DBHs
      ###__________________________________________________________###
      # Use the SfM-detected stems remaining after the filtering workflow
      # for the local DBH to height allometric relationship model.
        if(nrow(dbh_training_data_temp)>10){
          if(stringr::str_squish( tolower(dbh_model_local) ) == "rf"){
            # set random seed
            set.seed(21)

            ### tuning RF model
              # predictors and response to pass to randomForest functions
              predictors <- dbh_training_data_temp %>% dplyr::select(-c(treeID,stem_dbh_cm))
              response <- dbh_training_data_temp$stem_dbh_cm

              # implements steps to mitigate very long run-times when tuning random forests models
              optimal_mtry <- rf_tune_subsample(
                predictors = predictors
                , response = response
              )

              ### Run a randomForest model to predict HMD using various crown predictors
              # quiet this
              quiet_rf <- purrr::quietly(randomForest::randomForest)
              # run it
              stem_prediction_model <- quiet_rf(
                y = response
                , x = predictors
                , mtry = optimal_mtry
                , na.action = na.omit
              )

              # just get the result
              stem_prediction_model <- stem_prediction_model$result

            # stem_prediction_model
            # str(stem_prediction_model)

            # # variable importance plot
            #   randomForest::varImpPlot(stem_prediction_model, main = "RF variable importance plot for DBH estimate")

            ## Estimated versus observed DBH
            # data.frame(
            #   dbh_training_data_temp
            #   , predicted = stem_prediction_model$predicted
            # ) %>%
            # ggplot() +
            #   geom_abline() +
            #   geom_point(mapping = aes(x = stem_dbh_cm, y = predicted)) +
            #   scale_x_continuous(limits = c(0,max(dbh_training_data_temp$stem_dbh_cm)*1.05)) +
            #   scale_y_continuous(limits = c(0,max(dbh_training_data_temp$stem_dbh_cm)*1.05)) +
            #   labs(
            #     x = "SfM DBH (cm)"
            #     , y = "Predicted DBH (cm) by RF"
            #   ) +
            #   theme_light()

          }else{
            # population model with no random effects (i.e. no group-level variation)
            # Gamma distribution for strictly positive response variable dbh
            stem_prediction_model <- brms::brm(
              formula = stem_dbh_cm ~ 1 + tree_height_m
              , data = dbh_training_data_temp %>%
                  dplyr::select(stem_dbh_cm, tree_height_m)
              , family = brms::brmsfamily("Gamma", link = "log")
              , prior = c(brms::prior(gamma(0.01, 0.01), class = shape))
              , iter = 4000, warmup = 2000, chains = 4
              , cores = lasR::half_cores()
              , file = paste0(normalizePath(outfolder), "/local_dbh_height_model")
              , file_refit = "always"
            )
          }
          ###___________________________________________________________________###
          ### Predict missing DBH values for the top down crowns with no DBH ###
          ###___________________________________________________________________###
            crowns_sf_predict_only_temp <- tree_list %>%
              sf::st_drop_geometry() %>%
              dplyr::anti_join(
                dbh_training_data_temp %>%
                  dplyr::select(treeID)
                , by = dplyr::join_by("treeID")
              ) %>%
              dplyr::select(
                dbh_training_data_temp %>% dplyr::select(-c(stem_dbh_cm)) %>% names()
              )
            # str(crowns_sf_predict_only_temp)

            # get predicted dbh
            predicted_dbh_cm_temp <- predict(
              stem_prediction_model
              , crowns_sf_predict_only_temp %>% dplyr::select(-treeID)
            ) %>%
            dplyr::as_tibble() %>%
            dplyr::pull(1)
            # summary(predicted_dbh_cm_temp)

            ## combine predicted data with training data for full data set for all tree crowns with a matched tree top
            # nrow(crowns_sf)
            tree_tops <- tree_tops %>%
              # join training data
              dplyr::left_join(
                dbh_training_data_temp %>%
                  dplyr::mutate(is_training_data = T) %>%
                  dplyr::select(treeID, is_training_data, stem_dbh_cm)
                , by = dplyr::join_by("treeID")
              ) %>%
              # join with predicted data estimates
              dplyr::left_join(
                crowns_sf_predict_only_temp %>%
                  dplyr::mutate(
                    predicted_dbh_cm = predicted_dbh_cm_temp
                  ) %>%
                  dplyr::select(treeID, predicted_dbh_cm)
                , by = dplyr::join_by("treeID")
              ) %>%
              # clean up data and calculate metrics from dbh
              dplyr::mutate(
                is_training_data = dplyr::coalesce(is_training_data,F)
                , dbh_cm = dplyr::coalesce(stem_dbh_cm, predicted_dbh_cm, fia_est_dbh_cm)
              )
      }else{ # if(nrow(dbh_training_data_temp)>10)
        message(paste0(
          "Insufficient data to estimate DBH using trees provided in `treels_dbh_locations`..."
          , "\nReturning tree list with DBH estimates using FIA data instead."
        ))
      }
    } # end if treels_dbh_locations has data

  # fill in missing data if didn't go in to predict local dbh model
  if(max(grepl("is_training_data", names(tree_tops)))==0){
    tree_tops <- tree_tops %>%
      dplyr::mutate(
        is_training_data = F
        , stem_dbh_cm = as.numeric(NA)
        , predicted_dbh_cm = as.numeric(NA)
      )
  }

  # prep final data
    tree_tops <- tree_tops %>%
      dplyr::mutate(
        dbh_m = dbh_cm/100
        , radius_m = dbh_m/2
        , basal_area_m2 = pi * (radius_m)^2
        , basal_area_ft2 = basal_area_m2 * 10.764
        , ptcld_extracted_dbh_cm = stem_dbh_cm
        , ptcld_predicted_dbh_cm = predicted_dbh_cm
      ) %>%
      dplyr::select(-c(stem_dbh_cm, predicted_dbh_cm))
  # return
  if(treemap_data_finder_ans$which_treemap==2016){
    warning(paste0(
      "Treemap 2022 data has not been downloaded to package contents. You are currently using Treemap 2016."
      , "\n ... Use `get_treemap(force = T)` to update data"
    ))
  }
  return(tree_tops)
}
