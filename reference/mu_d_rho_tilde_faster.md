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

- fixed_P_multipliers:

  Numeric vector of length K containing values in (0,1\] equal to 1 -
  sum(fixed relative abundances in row k of P_tilde)

## Value

A derivative d mu_i / d rho_tilde_k
