# Compute the mean for each taxa across population

This function estimate the mean microbiome by averaging abundances
across individuals for each taxa. Optionally, it use a Dirichlet
distribution to simulate inter-individual variability centered arround
the mean.

## Usage

``` r
compute_mean_microbiome(microbiome, dir = F, n_ind = NULL, ao, mix.params)
```

## Arguments

- microbiome:

  A `matrix` of microbiome abundances (taxa in rows, individuals in
  columns).

- dir:

  A `logical`; if `TRUE`, use a Dirichlet distribution to simulate
  inter-individual variability. Default is `FALSE`.

- n_ind:

  The number of individuals to simulate with the Dirichlet distribution
  (required if `dir = TRUE`).

- ao:

  A numeric scalar used as the concentration parameter for the Dirichlet
  distribution.

- mix.params:

  A numeric vector (`length = 2`) specifying weights between Dirichlet
  samples and the original mean. The weigths should sum to 1.

## Value

A `numeric`vector of mean microbiome values if `dir = FALSE`. A `matrix`
of simulated microbiomes if `dir = TRUE`

## See also

[`compute_current_microbiome()`](compute_current_microbiome.md)
