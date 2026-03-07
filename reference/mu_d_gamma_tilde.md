# Calculate derivative of mu_ij with respect to an entry of gamma_tilde

Calculate derivative of mu_ij with respect to an entry of gamma_tilde

## Usage

``` r
mu_d_gamma_tilde(
  i,
  j,
  k_tilde,
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

- k_tilde:

  The element of gamma_tilde with respect to which to take derivative

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

- alpha_tilde:

  A numeric vector containing starting values of length \\M\\. If used,
  `Z_tilde_list` must be provided.

- Z_tilde_list:

  A list of length \\M + 1\\ containing matrices
  \\\tilde{Z}\_1,\dots,\tilde{Z}\_{M + 1}\\ to be linearly combined to
  generate `Z_tilde`: \\\tilde{Z} = \tilde{Z}\_{(1)} + \sum\_{m = 1}^M
  \exp(\tilde{\alpha}\_m)\tilde{Z}\_{(m + 1)}\\. If used, `alpha_tilde`
  must be provided.

- P_tilde:

  The spurious source relative abundance matrix (K_tilde x J)

- gamma_tilde:

  Spurious read intensity parameter

## Value

A derivative d mu_ij / d gamma_tilde_k_tilde
