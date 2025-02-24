#' @title implements steps to mitigate very long run-times when tuning random forests models
#'
#' @description
#' `rf_tune_subsample()` implements steps to mitigate very long run-times when tuning random forests models.
#' [randomForest::tuneRF()] enables model tuning by searching for the optimal `mtry` parameter (the number of variables randomly sampled as candidates at each split) using a cross-validation approach.
#' However, computational cost increases significantly with the number of observations as [randomForest::tuneRF()] performs cross-validation internally for each `mtry` value it tries.
#' With 100,000+ observations, each of these cross-validation runs involves building and evaluating many random forest trees, making the process very time-consuming.
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
#' @keywords internal
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
