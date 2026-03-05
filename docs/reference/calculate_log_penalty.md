# Calculate barrier penalty to add to objective function inside barrier algorithm

Calculate barrier penalty to add to objective function inside barrier
algorithm

## Usage

``` r
calculate_log_penalty(varying_lr_df, fixed_df, barrier_t)
```

## Arguments

- varying_lr_df:

  A data frame containing values of parameters that are treated as
  unknown, with relative abundance parameters represented on the log
  ratio scale (i.e., as phi and phi_tilde)

- fixed_df:

  A data frame containing values of parameters that are treated as known

- barrier_t:

  The current value of t, the barrier penalty parameter

## Value

The calculated value of the barrier penalty
