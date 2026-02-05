#' @title implements steps to mitigate very long run-times when tuning random forests models
#'
#' @description
#' `rf_tune_subsample()` implements steps to mitigate very long run-times when tuning random forests models.
#' [randomForest::tuneRF()] enables model tuning by searching for the optimal `mtry` parameter (the number of variables randomly sampled as candidates at each split) using a cross-validation approach.
#' However, computational cost increases significantly with the number of observations as [randomForest::tuneRF()] performs cross-validation internally for each `mtry` value it tries.
#' With 100,000+ observations, each of these cross-validation runs involves building and evaluating many random forest trees, making the process very time-consuming.
#'
#' The computational cost of random forests is driven by the repeated tree building process,
#' which involves recursive partitioning, bootstrapping, and feature subset selection.
#' These operations, when performed on massive datasets, result in a significant computational burden.
#'
#' `rf_tune_subsample()` remedies these issues via:
#'
#' * Reducing the `ntreeTry` parameter to a smaller value. Tuning will be less precise, but it will finish in a reasonable time. The `ntree` parameter can then be increased for the final model.
#' * Subsampling. Uses a smaller, representative subsample of the data (e.g., 10-20% of your data) to find a good `mtry` value on the subsample.
#'
#' @param predictors data.frame. predictor variable (x) data
#' @param response numeric. vector of response variable (y) data. observations should be ordered to match those in `predictors`
#' @param threshold numeric. the threshold number of observations, if observations exceed this threshold, subsampling is implemented
#' @param n_subsamples numeric. number of times to subsample and tune using [randomForest::tuneRF()]. The most common optimal `mtry` is returned from these subsample iterations.
#' @param ntree_try numeric. see [randomForest::tuneRF()]
#' @param step_factor numeric. see [randomForest::tuneRF()]
#' @param improve numeric. see [randomForest::tuneRF()]
#'
#' @return A numeric value to use in the `mtry` parameter of [randomForest::randomForest()]
#'
#' @noRd
#'
rf_tune_subsample <- function(
  predictors
  , response
  , threshold = 14444
  , n_subsamples = 4
  , ntree_try = 44
  , step_factor = 1
  , improve = 0.03
) {
  # set upper limit for ntreeTry
  max_ntreeTry <- 122
  # count the obs
  n_obs <- nrow(predictors)
  # set up list to store values from subsample runs
  optimal_mtry_values <- numeric(as.numeric(n_subsamples))
  # quiet the thing
  quiet_tuneRF <- purrr::quietly(randomForest::tuneRF)
  # if we don't need to subsample
  if(
    n_obs<=threshold
  ){
    # Tune on the full data if below threshold
    tune_result <- quiet_tuneRF(
      # randomForest::tuneRF(
      y = response
      , x = predictors
      , ntreeTry = max_ntreeTry # this is high for data with less obs
      , stepFactor = step_factor
      , improve = improve
      , plot = F
      , trace = F
    )
    # just get the result
    tune_result <- tune_result$result
    # Extract the optimal mtry value
    optimal_mtry <- tune_result %>%
      dplyr::as_tibble() %>%
      dplyr::filter(OOBError==min(OOBError)) %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::pull(mtry) %>%
      as.numeric()
  }else{
    # number for subsample
    sample_obs <- as.numeric(threshold)
    # run tuning on subsamples
    for (i in 1:n_subsamples) {
      # subsample
      subsample_indices <- sample(x = 1:n_obs, size = sample_obs, replace = F)
      subsample_predictors <- predictors[subsample_indices, ]
      subsample_response <- response[subsample_indices]
      # tune on subsample
      tune_result <- quiet_tuneRF(
        # randomForest::tuneRF(
        y = subsample_response
        , x = subsample_predictors
        , ntreeTry = min(as.numeric(ntree_try),max_ntreeTry)
        , stepFactor = step_factor
        , improve = improve
        , plot = F
        , trace = F
      )
      # just get the result
      tune_result <- tune_result$result
      # Extract the optimal mtry value
      optimal_mtry_values[i] <- tune_result %>%
        dplyr::as_tibble() %>%
        dplyr::filter(OOBError==min(OOBError)) %>%
        dplyr::filter(dplyr::row_number() == 1) %>%
        dplyr::pull(mtry) %>%
        as.numeric()
    }
    # Find most frequent optimal mtry
    optimal_mtry <- table(optimal_mtry_values) %>%
      sort(decreasing = TRUE) %>%
      .[1] %>%
      names() %>%
      as.numeric()
  }
  # ensure that the mtry value is not greater than the number of predictors
  optimal_mtry <- min(optimal_mtry, ncol(predictors))
  # return
  return(optimal_mtry)
}

################################################################################################################################################
# function to combine tuning and modelling
################################################################################################################################################
rf_tune_model <- function(
  predictors
  , response
  , tune_threshold = 14444
  , tune_n_subsamples = 4
  , tune_ntree_try = 44
  , tune_step_factor = 1
  , tune_improve = 0.03
  , ntree = 500
) {
  ### tuning RF model
  # implements steps to mitigate very long run-times when tuning random forests models
  optimal_mtry <- rf_tune_subsample(
    predictors = predictors
    , response = response
    , threshold = tune_threshold
    , n_subsamples = tune_n_subsamples
    , ntree_try = tune_ntree_try
    , step_factor = tune_step_factor
    , improve = tune_improve
  )

  ### Run a randomForest model
  # quiet this
  quiet_rf <- purrr::quietly(randomForest::randomForest)
  # run it
  mod <- quiet_rf(
    y = response
    , x = predictors
    , mtry = optimal_mtry
    , na.action = na.omit
    , ntree = ntree
  )

  # just get the result
  mod <- mod$result
  return(mod)
}
################################################################################################################################################
# function to subsample, tune, and model
################################################################################################################################################
rf_subsample_and_model <- function(
  predictors
  , response
  , tune_threshold = 14444
  , tune_n_subsamples = 4
  , tune_ntree_try = 44
  , tune_step_factor = 1
  , tune_improve = 0.03
  , mod_n_subsample = 20000
  , ntree = 500
) {
  if(length(response)<=mod_n_subsample){
    mod <- rf_tune_model(
        predictors = predictors
        , response = response
        , tune_threshold = tune_threshold
        , tune_n_subsamples = tune_n_subsamples
        , tune_ntree_try = tune_ntree_try
        , tune_step_factor = tune_step_factor
        , tune_improve = tune_improve
        , ntree = ntree
      )
  }else{
    subsample_indices <- sample(x = 1:length(response), size = mod_n_subsample, replace = F)
    subsample_y <- response[subsample_indices]
    subsample_x <- predictors[subsample_indices,]
    mod <- rf_tune_model(
        predictors = subsample_x
        , response = subsample_y
        , tune_threshold = tune_threshold
        , tune_n_subsamples = tune_n_subsamples
        , tune_ntree_try = tune_ntree_try
        , tune_step_factor = tune_step_factor
        , tune_improve = tune_improve
        , ntree = ntree
      )
  }
  return(mod)
}
################################################################################################################################################
# function to subsample, tune, and model N TIMES
# iterates the subsampling and model fitting process n_iterations times
# stores the predictions from each model in a list
################################################################################################################################################
rf_subsample_and_model_n_times <- function(
  predictors
  , response
  , tune_threshold = 14444
  , tune_n_subsamples = 4
  , tune_ntree_try = 44
  , tune_step_factor = 1
  , tune_improve = 0.03
  , mod_n_subsample = 20000
  , mod_n_times = 3
  , ntree = 500
) {
  # set up blank list
  mod_list <- list()
  if(length(response)<=mod_n_subsample){
    mod_list[[1]] <- rf_subsample_and_model(
        predictors = predictors
        , response = response
        , tune_threshold = tune_threshold
        , tune_n_subsamples = tune_n_subsamples
        , tune_ntree_try = tune_ntree_try
        , tune_step_factor = tune_step_factor
        , tune_improve = tune_improve
        , mod_n_subsample = mod_n_subsample
        , ntree = ntree
      )
  }else{
    for (i in 1:mod_n_times) {
      mod_list[[i]] <- rf_subsample_and_model(
          predictors = predictors
          , response = response
          , tune_threshold = tune_threshold
          , tune_n_subsamples = tune_n_subsamples
          , tune_ntree_try = tune_ntree_try
          , tune_step_factor = tune_step_factor
          , tune_improve = tune_improve
          , mod_n_subsample = mod_n_subsample
          , ntree = ntree
        )
    }
  }
  return(mod_list)
}
################################################################################################################################################
# function to model average predictions using model list
# Model averaging can improve the robustness and accuracy of random forest models, especially when dealing with large datasets
# predict_df must have predictor (x) columns used to develop model
################################################################################################################################################
rf_model_avg_predictions <- function(mod_list, predict_df){
  # get model predictions for each model
  predicted_temp <-
    1:length(mod_list) %>%
      purrr::map(function(x){
        predict(
          mod_list[[x]]
          , predict_df
        ) %>%
        dplyr::as_tibble() %>%
        dplyr::pull(1) %>%
        dplyr::as_tibble() %>%
        dplyr::rename(value=1) %>%
        dplyr::mutate(id=dplyr::row_number())
      }) %>%
      dplyr::bind_rows()
  # model average
  predicted_temp <- predicted_temp %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(
      predicted = mean(value, na.rm = T)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(id) %>%
    dplyr::select(predicted)
  # return
  return(predicted_temp)
}
