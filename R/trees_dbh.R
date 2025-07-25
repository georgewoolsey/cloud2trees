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
#' @param dbh_model string. Set the model to use for local dbh-height allometry. Can be "rf" for random forest or "lin" for linear
#' @param treels_dbh_locations sf. Return from [treels_stem_dbh()].
#' Must also provide crown polygons (as returned from [raster2trees()]) in the `tree_list` data as an `sf` class object with POLYGON geometry (see [sf::st_geometry_type()])
#' If a valid file is provided, will make DBH predictions based on this training data instead of from the regional model from the FIA data
#' @param boundary_buffer numeric. Set the buffer (m) for the study area boundary to filter the FIA plot data based on TreeMap
#' @param input_treemap_dir directory where Treemap 2016 exists. Use [get_treemap()] first.
#' @param outfolder string. The path of a folder to write the model data to
#'
#' @references
#' * [https://doi.org/10.2737/RDS-2021-0074](https://doi.org/10.2737/RDS-2021-0074)
#' Riley, Karin L.; Grenfell, Isaac C.; Finney, Mark A.; Shaw, John D. 2021. TreeMap 2016: A tree-level model of the forests of the conterminous United States circa 2016. Fort Collins, CO: Forest Service Research Data Archive.
#'
#' * [https://doi.org/10.3390/f13122077](https://doi.org/10.3390/f13122077)
#' Tinkham et al. (2022). Modeling the missing DBHs: Influence of model form on UAV DBH characterization. Forests, 13(12), 2077.
#'
#' @return Returns a spatial data frame of individual trees.
#'
#' @examples
#'  \dontrun{
#'  # example tree list
#'  tl <- dplyr::tibble(
#'      treeID = c(1:21)
#'      , tree_x = rnorm(n=21, mean = 458064, sd = 11)
#'      , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
#'      , tree_height_m = exp(rgamma(n = 21, shape = (7/4)^2, rate = (4^2)/7))
#'    )
#'  # call the function
#'  tl_dbh <- trees_dbh(tree_list = tl, crs = "32613")
#'  # what?
#'  tl_dbh %>% class()
#'  tl_dbh %>% dplyr::select(tidyselect::contains("dbh_cm")) %>% dplyr::glimpse()
#'  tl_dbh %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=dbh_cm))
#'  }
#' @export
#'
trees_dbh <- function(
  tree_list
  , crs = NA
  , study_boundary = NA
  , dbh_model = "lin"
  , treels_dbh_locations = NA
  , boundary_buffer = 50
  , input_treemap_dir = NULL
  , outfolder = tempdir()
) {
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
      )

    # ggplot(treemap_rast %>% as.data.frame(xy=T) %>% rename(f=3)) +
    #   geom_tile(aes(x=x,y=y,fill=as.factor(f))) +
    #   scale_fill_viridis_d(option = "turbo") +
    #   theme_light() + theme(legend.position = "none")

    ### get weights for weighting each tree in the population models
    # treemap id = tm_id for linking to tabular data
    tm_id_weight_temp <- terra::freq(treemap_rast) %>%
      dplyr::select(-layer) %>%
      dplyr::rename(tm_id = value, tree_weight = count) %>%
      dplyr::mutate(tm_id = as_character_safe(tm_id))
    # str(tm_id_weight_temp)

    ### get the TreeMap FIA tree list for only the plots included
    treemap_trees_df <- readr::read_csv(
        file.path(find_ext_data_ans$treemap_dir, "treemap2016_tree_table.csv")
        , col_select = c(
          tm_id
          , CN
          , SPECIES_SYMBOL
          , STATUSCD
          , DIA
          , HT
        )
        , progress = F
        , show_col_types = F
      ) %>%
      dplyr::rename_with(tolower) %>%
      dplyr::mutate(
        cn = as_character_safe(cn)
        , tm_id = as_character_safe(tm_id)
      ) %>%
      dplyr::left_join(
        tm_id_weight_temp
        , by = dplyr::join_by("tm_id")
      ) %>%
      dplyr::left_join(
        tm_id_weight_temp %>% dplyr::rename(cn = tm_id)
        , by = dplyr::join_by("cn")
      ) %>%
      dplyr::mutate(tree_weight = dplyr::coalesce(tree_weight.x, tree_weight.y)) %>%
      dplyr::select(-c(tree_weight.x, tree_weight.y)) %>%
      dplyr::filter(
        # keep live trees only: 1=live;2=dead
        statuscd == 1
        & !is.na(dia)
        & !is.na(ht)
        & !is.na(tree_weight)
      ) %>%
      dplyr::mutate(
        dbh_cm = dia*2.54
        , tree_height_m = ht/3.28084
      ) %>%
      dplyr::select(-c(statuscd,dia,ht))

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
          treemap_trees_df
          , paste0(normalizePath(outfolder), "/regional_dbh_height_model_training_data.csv")
          , row.names = F
        )

    ###__________________________________________________________###
    ### Regional model of DBH as predicted by height
    ### population model of dbh on height, non-linear
    ### used to filter sfm dbhs
    ###__________________________________________________________###
    # population model with no random effects (i.e. no group-level variation)
    # non-linear model form with Gamma distribution for strictly positive response variable dbh
    # set up prior
    p_temp <- brms::prior(normal(1, 2), nlpar = "b1") +
      brms::prior(normal(0, 2), nlpar = "b2")
    mod_nl_pop <- brms::brm(
      formula = brms::bf(
        formula = dbh_cm|weights(tree_weight) ~ (b1 * tree_height_m) + tree_height_m^b2
        , b1 + b2 ~ 1
        , nl = TRUE # !! specify non-linear
      )
      , data = treemap_trees_df
      , prior = p_temp
      , family = brms::brmsfamily("Gamma")
      , iter = 4000, warmup = 2000, chains = 4
      , cores = lasR::half_cores()
      , file = paste0(normalizePath(outfolder), "/regional_dbh_height_model")
      , file_refit = "always"
    )
    # plot(mod_nl_pop)
    # summary(mod_nl_pop)

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
          if(tolower(dbh_model) == "rf"){
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
  return(tree_tops)
}
