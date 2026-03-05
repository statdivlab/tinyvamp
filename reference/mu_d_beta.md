# Calculate derivative of mu_ij with respect to B

Calculate derivative of mu_ij with respect to B

## Usage

``` r
mu_d_beta(
  i,
  j,
  q,
  gammas,
  B,
  X,
  Z,
  P,
  X_tilde,
  Z_tilde,
  Z_tilde_gamma_cols,
  alpha_tilde = NULL,
  Z_tilde_list = NULL,
  P_tilde,
  gamma_tilde
)
```

## Arguments

- i:

  The sample index (must be in 1, ..., n)

- j:

  The taxon index (must be in 1, ..., J)

- q:

  Which row of B to take derivative with respect to (must be in 1, ...,
  p)

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

A derivative d mu_ij / d B_qj
