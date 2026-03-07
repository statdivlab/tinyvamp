# Apply the Bayesian subsampled bootstrap to a fitted tinyvamp model

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

## Arguments

- W:

  A numeric matrix of observed counts with \\n\\ rows (samples) and
  \\J\\ columns (taxa).

- fitted_model:

  A fitted tinyvamp model object, typically returned by
  [`estimate_parameters()`](https://statdivlab.github.io/tinyvamp/reference/estimate_parameters.md),
  containing the original fit and all inputs needed for refitting during
  bootstrap iterations.

- n_boot:

  Integer; the number of bootstrap replicates.

- m:

  Optional numeric subsample size parameter for the Bayesian subsampled
  bootstrap. If `NULL`, defaults to `sqrt(n)`, where \\n\\ is the number
  of rows of `W`.

- alpha:

  Significance level used to construct two-sided confidence intervals.
  The default `0.05` yields 95\\

  parallelizeLogical; if `TRUE`, bootstrap replicates are fit in
  parallel using
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html). If
  `FALSE`, replicates are fit sequentially.

  ncoresInteger; number of cores to use when `parallelize = TRUE`.

  seedOptional integer random seed. If `NULL`, a seed value of `0` is
  used.

  return_modelsLogical; if `TRUE`, return the fitted bootstrap model
  objects in addition to the confidence interval summary.

  verboseLogical; if `TRUE`, print progress messages during sequential
  bootstrap fitting and pass verbose output to the internal model
  fitting routine.

  adjustLogical; if `TRUE`, apply a finite-sample adjustment to the
  bootstrap deviations based on the number of non-`"gamma"` varying
  parameters.

A list with components:

- `ci`:

  A data frame based on `fitted_model$varying` with added columns
  `lower_ci` and `upper_ci`.

- `bootstrapped_models`:

  If `return_models = TRUE`, a list of fitted model objects from the
  bootstrap replicates; otherwise `NULL`.

Compute bootstrap confidence intervals for the varying parameters of a
fitted tinyvamp model using a Bayesian subsampled bootstrap procedure.
Each bootstrap replicate refits the model using randomly generated
observation weights, and confidence intervals are formed from the
empirical quantiles of the bootstrap distribution. For each bootstrap
replicate, the function generates a vector of Bayesian bootstrap weights
from a gamma distribution, rescales them, and uses them as observation
weights in a call to
[`estimate_parameters()`](https://statdivlab.github.io/tinyvamp/reference/estimate_parameters.md).
The resulting bootstrap distribution is centered at the original fitted
values and scaled by \\\sqrt{m}\\. Confidence intervals are then
obtained by inverting the bootstrap quantiles and scaling by \\1 /
\sqrt{n}\\.When `adjust = TRUE`, the bootstrap deviations are multiplied
by a correction factor \$\$\sqrt{\frac{nJ}{nJ - p}}\$\$ where \\p\\ is
the number of varying parameters whose `param` field is not equal to
`"gamma"`.Confidence interval bounds for parameters labeled `"P"` or
`"P_tilde"` are truncated to lie in \\\[0, 1\]\\.
