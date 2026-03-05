# Fit tinyvamp model to HTS microbiome data

This function fits a model to HTS microbiome data that allows for
estimation of detection efficiency effects as well as modeling of
spurious read sources (e.g., contamination).

## Usage

``` r
estimate_parameters(
  W,
  X,
  Z,
  Z_tilde = NULL,
  Z_tilde_gamma_cols,
  gammas,
  gammas_fixed_indices,
  P,
  P_fixed_indices,
  B,
  B_fixed_indices,
  X_tilde,
  P_tilde,
  P_tilde_fixed_indices,
  gamma_tilde,
  gamma_tilde_fixed_indices,
  alpha_tilde = NULL,
  Z_tilde_list = NULL,
  barrier_t = 1,
  barrier_scale = 10,
  max_barrier = 1e+12,
  initial_conv_tol = 1000,
  final_conv_tol = 0.1,
  constraint_tolerance = 1e-10,
  hessian_regularization = 0.01,
  criterion = "Poisson",
  profile_P = TRUE,
  barrier_maxit = 500,
  profiling_maxit = 25,
  wts = NULL,
  verbose = FALSE,
  bootstrap_failure_cutoff = NULL,
  tinker_zeroes = 0.1,
  return_variance = FALSE
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

- gammas:

  A numeric vector of length n of starting values for read intensity
  parameter gamma

- gammas_fixed_indices:

  A logical vector of length n whose \\i\\-th entry is TRUE if the
  \\i\\-th entry of gamma should be treated as fixed and known, and
  FALSE otherwise

- P:

  A \\K \times J\\ numeric matrix giving initial values for the relative
  abundance matrix.

- P_fixed_indices:

  P_fixed_indices A \\K \times J\\ logical matrix specifying any entries
  of P that are known. If known, the corresponding values from `P` will
  be treated as the fixed, known values.

- B:

  A \\p \times J\\ numeric matrix giving initial values for the sample
  efficiencies.

- B_fixed_indices:

  A \\p \times J\\ logical matrix specifying any entries of B that are
  known. If known, the corresponding values from `B` will be treated as
  the fixed, known values.

- X_tilde:

  A \\\tilde{K} \times p\\ matrix giving the spurious read source
  efficiency design matrix

- P_tilde:

  A \\\tilde{K} \times J\\ numeric matrix giving initial values for the
  spurious read source relative abundances.

- P_tilde_fixed_indices:

  A \\\tilde{K} \times J\\ logical matrix indicating if the \\(i,j)\\th
  entry of `P_tilde` should be treated as fixed and known.

- gamma_tilde:

  A numeric vector of length \\\tilde{K}\\ of starting values for
  spurious read intensity parameter gamma_tilde

- gamma_tilde_fixed_indices:

  A logical vector of length \\\tilde{K}\\ whose \\i\\-th entry is TRUE
  if the \\i\\-th entry of gamma_tilde should be treated as fixed and
  known, and FALSE otherwise.

- alpha_tilde:

  A numeric vector containing starting values of length \\M\\. If used,
  `Z_tilde_list` must be provided.

- Z_tilde_list:

  A list of length \\M + 1\\ containing matrices
  \\\tilde{Z}\_1,\dots,\tilde{Z}\_{M + 1}\\ to be linearly combined to
  generate `Z_tilde`: \\\tilde{Z} = \tilde{Z}\_{(1)} + \sum\_{m = 1}^M
  \exp(\tilde{\alpha}\_m)\tilde{Z}\_{(m + 1)}\\. If used, `alpha_tilde`
  must be provided.

- barrier_t:

  Starting value of reciprocal barrier penalty coef. Defaults to 1.

- barrier_scale:

  Increments for value of barrier penalty. Defaults to 10.

- max_barrier:

  Maximum value of barrier_t. Defaults to 1e12.

- constraint_tolerance:

  The tolerance for the augmented Lagrangian algorithm. Final estimates
  of P are relative abundances to within `constraint_tolerance` of 1,
  i.e., abs(sum p_kj - 1) \< `constraint_tolerance`. Defaults to 1e-10.

- hessian_regularization:

  The second step of optimization involves a quadratic approximation to
  the likelihood, for which we use a modified Taylor series for
  stability. This is the constant that dampens the second term. Defaults
  to 0.01.

- criterion:

  Should the algorithm return the Poisson maximum likelihood estimates
  or the reweighted Poisson maximum likelihood estimates? Options are
  "Poisson" or "reweighted_Poisson".

- profile_P:

  Defaults to TRUE Run profiling step after barrier algorithm has run?
  If TRUE, this step is performed, possibly setting some estimated
  relative abundances in P equal to zero. If FALSE, profiling step is
  skipped and back-transformed log-ratio parameter estimated via barrier
  algorithm is returned for P.

- barrier_maxit:

  The maximum number of iterations for the barrier method

- profiling_maxit:

  Maximum number of iterations to run profiling step in P for (default
  is 25).

- wts:

  Weights for reweighting the likelihood contributions. This is usually
  done to improve efficiency. Defaults to NULL. We compute the weights
  for you even if you choose `criterion = "reweighted_Poisson"`.

- verbose:

  Do you want to know what I'm doing? Defaults to FALSE.

- bootstrap_failure_cutoff:

  Defaults to NULL.

- tinker_zeroes:

  Because the barrier method can only be applied to relative abundances
  in the interior of the simplex, tinker_zeroes divided by the number of
  taxa is added to all relative abundances to be estimated before the
  barrier method is applied. Default 0.1.

- return_variance:

  Defaults to FALSE.

## Value

A list containing estimated parameter values, along with the given
inputs

## Author

David Clausen

Amy Willis
