# implements steps to mitigate very long run-times when tuning random forests models

`rf_tune_subsample()` implements steps to mitigate very long run-times
when tuning random forests models.
[`randomForest::tuneRF()`](https://rdrr.io/pkg/randomForest/man/tuneRF.html)
enables model tuning by searching for the optimal `mtry` parameter (the
number of variables randomly sampled as candidates at each split) using
a cross-validation approach. However, computational cost increases
significantly with the number of observations as
[`randomForest::tuneRF()`](https://rdrr.io/pkg/randomForest/man/tuneRF.html)
performs cross-validation internally for each `mtry` value it tries.
With 100,000+ observations, each of these cross-validation runs involves
building and evaluating many random forest trees, making the process
very time-consuming.

The computational cost of random forests is driven by the repeated tree
building process, which involves recursive partitioning, bootstrapping,
and feature subset selection. These operations, when performed on
massive datasets, result in a significant computational burden.

`rf_tune_subsample()` remedies these issues via:

- Reducing the `ntreeTry` parameter to a smaller value. Tuning will be
  less precise, but it will finish in a reasonable time. The `ntree`
  parameter can then be increased for the final model.

- Subsampling. Uses a smaller, representative subsample of the data
  (e.g., 10-20% of your data) to find a good `mtry` value on the
  subsample.

## Usage

``` r
rf_tune_subsample(
  predictors,
  response,
  threshold = 14444,
  n_subsamples = 4,
  ntree_try = 44,
  step_factor = 1,
  improve = 0.03
)
```

## Arguments

- predictors:

  data.frame. predictor variable (x) data

- response:

  numeric. vector of response variable (y) data. observations should be
  ordered to match those in `predictors`

- threshold:

  numeric. the threshold number of observations, if observations exceed
  this threshold, subsampling is implemented

- n_subsamples:

  numeric. number of times to subsample and tune using
  [`randomForest::tuneRF()`](https://rdrr.io/pkg/randomForest/man/tuneRF.html).
  The most common optimal `mtry` is returned from these subsample
  iterations.

- ntree_try:

  numeric. see
  [`randomForest::tuneRF()`](https://rdrr.io/pkg/randomForest/man/tuneRF.html)

- step_factor:

  numeric. see
  [`randomForest::tuneRF()`](https://rdrr.io/pkg/randomForest/man/tuneRF.html)

- improve:

  numeric. see
  [`randomForest::tuneRF()`](https://rdrr.io/pkg/randomForest/man/tuneRF.html)

## Value

A numeric value to use in the `mtry` parameter of
[`randomForest::randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html)
