# Convert parameter values stored in data frame format to matrix format

Convert parameter values stored in data frame format to matrix format

## Usage

``` r
dataframes_to_parameters(fixed_df, varying_df)
```

## Arguments

- fixed_df:

  A dataframe containing values of model parameters treated as fixed and
  known (i.e. held constant at known values)

- varying_df:

  A dataframe containing current values of model parameters treated as
  fixed and unknown (i.e., parameters to be estimated)

## Value

A list containing

- P:

  Specimen relative abundance matrix (of dimension K x J)

- P_tilde:

  Spurious read source relative abundance matrix (of dimension K-tilde x
  J)

- B:

  A matrix of detection efficiencies (of dimension p x J)

- gammas:

  An n-vector of sample-specific read intensities

- gamma_tilde:

  A-vector of spurious read source intensities (of length K-tilde)

## Author

David Clausen
