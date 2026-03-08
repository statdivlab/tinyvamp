# Calculate derivative of mu_ij with respect to an entry of P

Calculate derivative of mu_ij with respect to an entry of P

## Usage

``` r
mu_d_P(i, j, m, gammas, B, X, Z, P)
```

## Arguments

- i:

  The sample index (must be in 1, ..., n)

- j:

  The taxon index (must be in 1, ..., J)

- m:

  The row of P with respect to which to take derivative

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

## Value

A derivative d mu_ij / d P_kj
