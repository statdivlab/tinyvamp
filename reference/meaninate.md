# Calculate the mean given parameters

Calculate the mean given parameters

## Usage

``` r
meaninate(
  gammas,
  B,
  X,
  Z,
  P,
  X_tilde,
  Z_tilde = NULL,
  Z_tilde_gamma_cols,
  P_tilde,
  gamma_tilde,
  alpha_tilde = NULL,
  Z_tilde_list = NULL,
  return_separate = FALSE,
  exclude_gammas = FALSE
)
```

## Arguments

- gammas:

  A numeric vector of length n of starting values for read intensity
  parameter gamma

- B:

  A \\p \times J\\ numeric matrix giving initial values for the sample
  efficiencies.

- X:

  The sample efficiency design – an \\n \times p\\ matrix

- Z:

  The sample-specimen design – an \\n \times K\\ matrix whose \\ij\\-th
  entry indicates the proportional contribution of specimen \\j\\ to
  sample \\i\\. Rows must sum to 1 or be identically 0.

- P:

  A \\K \times J\\ numeric matrix giving initial values for the relative
  abundance matrix.

- X_tilde:

  A \\\tilde{K} \times p\\ matrix giving the spurious read source
  efficiency design matrix

- Z_tilde:

  The spurious read design – an \\n \times \tilde{K}\\ matrix where
  \\\tilde{K}\\ is the number of spurious read sources modeled.

- Z_tilde_gamma_cols:

  A numeric vector containing the columns of Z_tilde which should be
  multiplied by exp(gamma).

- P_tilde:

  A \\\tilde{K} \times J\\ numeric matrix giving initial values for the
  spurious read source relative abundances.

- gamma_tilde:

  A numeric vector of length \\\tilde{K}\\ of starting values for
  spurious read intensity parameter gamma_tilde

- alpha_tilde:

  A numeric vector containing starting values of length \\M\\. If used,
  `Z_tilde_list` must be provided.

- Z_tilde_list:

  A list of length \\M + 1\\ containing matrices
  \\\tilde{Z}\_1,\dots,\tilde{Z}\_{M + 1}\\ to be linearly combined to
  generate `Z_tilde`: \\\tilde{Z} = \tilde{Z}\_{(1)} + \sum\_{m = 1}^M
  \exp(\tilde{\alpha}\_m)\tilde{Z}\_{(m + 1)}\\. If used, `alpha_tilde`
  must be provided.

- return_separate:

  Boolean. Return the summed mean, or separate the sample and
  contamination pieces. Defaults to FALSE.

- exclude_gammas:

  Boolean, defaults to FALSE. Should the gamma components be ignored?
