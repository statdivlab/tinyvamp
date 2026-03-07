# Criterion Evaluation Function

Evaluates the objective function used during model fitting under either
a Poisson likelihood criterion or a GMM criterion.

## Usage

``` r
evaluate_criterion_lr(
  W,
  X,
  Z,
  Z_tilde,
  Z_tilde_gamma_cols,
  Z_tilde_list = NULL,
  X_tilde,
  fixed_df,
  varying_df,
  varying_lr_df = NULL,
  barrier_t = NULL,
  criterion = "Poisson",
  lr_scale = TRUE,
  include_log_penalty = TRUE,
  wts = NULL,
  gmm_inv_wts = NULL,
  return_gmm_inv_weights = FALSE
)
```

## Arguments

- W:

  An \\n \times J\\ matrix of numeric HTS output (e.g., read counts,
  coverages, etc.)

- X:

  The sample efficiency design – an \\n \times p\\ matrix

- Z:

  The sample-specimen design – an \\n \times K\\ matrix whose \\ij\\-th
  entry indicates the proportional contribution of specimen \\j\\ to
  sample \\i\\. Rows must sum to 1 or be identically 0.

- Z_tilde:

  The spurious read design – an \\n \times \tilde{K}\\ matrix where
  \\\tilde{K}\\ is the number of spurious read sources modeled.

- Z_tilde_gamma_cols:

  A numeric vector containing the columns of Z_tilde which should be
  multiplied by exp(gamma).

- Z_tilde_list:

  Optional list-form representation of `Z_tilde` or related auxiliary
  structures used by `meaninate()`.

- X_tilde:

  Design matrix associated with spurious-read or auxiliary components of
  the mean model.

- fixed_df:

  A data frame of fixed model parameters.

- varying_df:

  A data frame of parameters currently being optimized, represented on
  the natural scale unless `lr_scale = TRUE`.

- varying_lr_df:

  Optional data frame of varying parameters on the log-ratio scale. Used
  when `lr_scale = TRUE`.

- barrier_t:

  Optional numeric tuning parameter for the log-ratio barrier penalty.

- criterion:

  Character string specifying the criterion to evaluate. Currently one
  of `"Poisson"` or `"GMM"`.

- lr_scale:

  Logical; if `TRUE`, convert log-ratio scale parameters to relative
  abundance scale before evaluation.

- include_log_penalty:

  Logical; if `TRUE`, include the log-barrier penalty when
  `lr_scale = TRUE`.

- wts:

  Optional numeric weights passed to `poisson_criterion()`.

- gmm_inv_wts:

  Optional inverse weighting matrix or vector for the GMM criterion. If
  `NULL`, it is estimated internally using `get_gmm_inv_weights()`.

- return_gmm_inv_weights:

  Logical; if `TRUE` and `criterion = "GMM"`, return both the GMM
  criterion value and the inverse weights.

## Value

If `criterion = "Poisson"`, a single numeric criterion value.

If `criterion = "GMM"` and `return_gmm_inv_weights = FALSE`, a single
numeric criterion value.

If `criterion = "GMM"` and `return_gmm_inv_weights = TRUE`, a list with
components:

- gmm_crit:

  The numeric GMM criterion value.

- inv_wts:

  The inverse weights used in the GMM calculation.
