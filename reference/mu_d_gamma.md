# Calculate derivative of mu_ij with respect to ith entry of gamma

Calculate derivative of mu_ij with respect to ith entry of gamma

## Usage

``` r
mu_d_gamma(
  i,
  j,
  gammas,
  B,
  X,
  Z,
  P,
  X_tilde,
  Z_tilde,
  Z_tilde_gamma_cols,
  P_tilde,
  gamma_tilde
)
```

## Arguments

- i:

  The sample index (must be in 1, ..., n)

- j:

  The taxon index (must be in 1, ..., J)

- gammas:

  Numeric vector of read intensities

- B:

  Detection efficiency matrix

- X:

  The efficiency design matrix (n x p)

- Z:

  The sample design matrix (n x K)

- P:

  The sample relative abundance matrix (K x J)

- X_tilde:

  The spurious read efficiency design (K_tilde x p)

- Z_tilde:

  The spurious read design (n x K_tilde)

- Z_tilde_gamma_cols:

  Numeric vector containing indexes of columns of Z_tilde to scale by
  exp(gamma); NULL if no columns to be scaled

- P_tilde:

  The spurious source relative abundance matrix (K_tilde x J)

- gamma_tilde:

  Spurious read intensity parameter

## Value

A derivative d mu_ij / d gamma_i
