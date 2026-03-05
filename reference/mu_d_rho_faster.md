# Calculate derivative of mu_ij with respect to a row of matrix-valued parameter rho

Calculate derivative of mu_ij with respect to a row of matrix-valued
parameter rho

## Usage

``` r
mu_d_rho_faster(
  i,
  J,
  k,
  gammas,
  B,
  X,
  Z,
  Ak_list,
  rho_k,
  fixed_P_multipliers,
  proportion_scale = FALSE
)
```

## Arguments

- i:

  The sample index (must be in 1, ..., n)

- J:

  The total number of taxa modeled

- k:

  Row index (which row of rho with respect to which to take derivative)

- gammas:

  Numeric vector of read intensities

- B:

  Detection efficiency matrix

- X:

  The efficiency design matrix (n x p)

- Z:

  The sample design matrix (n x K)

- Ak_list:

  List containing matrices that map back-transformed rho to entries of P

- rho_k:

  Value of kth row of rho

- fixed_P_multipliers:

  Numeric vector of length K containing values in (0,1\] equal to 1 -
  sum(fixed relative abundances in row k of P) Z_tilde to scale by
  exp(gamma); NULL if no columns to be scaled

## Value

A derivative d mu_i / d rho_k
