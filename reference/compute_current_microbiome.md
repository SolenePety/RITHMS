# Compute microbiomes for all individuals of current generation gathering all relevant object already computed.

Compute microbiomes for all individuals of current generation gathering
all relevant object already computed.

## Usage

``` r
compute_current_microbiome(
  beta,
  current_genotypes,
  mother_microbiomes,
  mean_microbiome,
  noise = 0.1,
  taxa_scale,
  lambda = 0.5,
  dir = F,
  thetaX
)
```

## Arguments

- beta:

  A `matrix` of taxa (rows) by SNPs (columns) representing the
  multiplicative effect of SNPs on taxa abundances. Typically the output
  of [`compute_beta_matrix_cluster()`](compute_beta_matrix_cluster.md)

- current_genotypes:

  A `matrix` of genotypes for the individuals of the current generation.

- mother_microbiomes:

  A `matrix` containing the CLR-transformed abundances from the
  microbiome of the current mothers.

- mean_microbiome:

  A `numeric` vector representing the average abundances of the
  microbiome within the population.

- noise:

  A `numeric` scalar indicating the standard deviation of the Gaussian
  noise to add during beta g construction
  [`compute_beta_g()`](compute_beta_g.md)(multiplicative genetic effect
  on taxa abundances).

- taxa_scale:

  A `numeric` scalar to scale the noise added to each OTUs abundance,
  see [`compute_beta_g()`](compute_beta_g.md).

- lambda:

  A `numeric` between 0 and 1 controlling the relative contribution of
  the mother microbiome compared to the average microbiome.

- dir:

  A `logical`;

- thetaX:

  A `matrix` representing the environmental effect (optional).

## Value

A `matrix` of CLR-transformed abundances of the microbiome for each
individual of the current generation.

## See also

[`compute_mean_microbiome()`](compute_mean_microbiome.md)
