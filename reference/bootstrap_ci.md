# Apply the Bayesian subsampled bootstrap to a fitted tinyvamp model

Apply the Bayesian subsampled bootstrap to a fitted tinyvamp model

## Usage

``` r
bootstrap_ci(
  W,
  fitted_model,
  n_boot,
  m = NULL,
  alpha = 0.05,
  parallelize = FALSE,
  ncores = 5,
  seed = NULL,
  return_models = FALSE,
  verbose = FALSE,
  adjust = FALSE
)
```
