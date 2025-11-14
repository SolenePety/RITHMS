# Compute root value for beta matrix construction

This function computes the parameter used to build a beta matrix with a
given correlation structure between taxa. It is used in the beta matrix
construction.

## Usage

``` r
root(rho)
```

## Arguments

- rho:

  A numeric value between 0 and 1 representing the target correlation
  between taxa within the same cluster. This parameter is used to induce
  correlation in the beta matrix.
