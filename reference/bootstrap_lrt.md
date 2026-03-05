# Apply a likelihod ratio test

Apply a likelihod ratio test

## Usage

``` r
bootstrap_lrt(
  W,
  fitted_model,
  null_param,
  n_boot,
  m = NULL,
  recalculate_W0 = FALSE,
  parallelize = FALSE,
  ncores = 1,
  save_models = FALSE,
  verbose = FALSE
)
```

## Arguments

- W:

  An \\n \times J\\ matrix of numeric HTS output (e.g., read counts,
  coverages, etc.)

- fitted_model:

  The model fitted under the alternative

- null_param:

  The model fitted under the null. This could be a modified copy of the
  model fit under the alternative.

- n_boot:

  The number of bootstrap resamples to perform

- m:

  description

- recalculate_W0:

  Boolean. Should rescaling W. TODO. Defaults to FALSE.

- parallelize:

  Boolean. do you want to parallelize it?

- ncores:

  the number of cores to parallelize over

- save_models:

  TODO

- verbose:

  Do you want to know what I'm doing? Defaults to FALSE.
