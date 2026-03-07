# Calculate derivative of mu_i with respect to a row of matrix-valued parameter rho_tilde

Calculate derivative of mu_i with respect to a row of matrix-valued
parameter rho_tilde

## Usage

``` r
mu_d_rho_tilde_faster(
  i,
  J,
  k_tilde,
  gammas,
  B,
  rho_tilde_k,
  A_tilde_k_list,
  fixed_P_tilde_multipliers,
  X_tilde,
  Z_tilde,
  Z_tilde_gamma_cols,
  alpha_tilde = NULL,
  Z_tilde_list = NULL,
  gamma_tilde,
  proportion_scale = FALSE
)
```

## Arguments

- i:

  The sample index (must be in 1, ..., n)

- J:

  The total number of taxa modeled

- k_tilde:

  Row index (which row of rho_tilde with respect to which to take
  derivative)

- gammas:

  Numeric vector of read intensities

- B:

  Detection efficiency matrix

- rho_tilde_k:

  Value of kth row of rho_tilde

- A_tilde_k_list:

  List containing matrices that map back-transformed rho_tilde to
  entries of P_tilde

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

- fixed_P_multipliers:

  Numeric vector of length K containing values in (0,1\] equal to 1 -
  sum(fixed relative abundances in row k of P_tilde)

## Value

A derivative d mu_i / d rho_tilde_k
