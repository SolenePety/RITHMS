# Compute product of matrixes based on few parameters

This function compute the product between the matrix of OTUs-specific
genetic effects and the matrix of genotypes. This is part of
[`gen_effect_calibration()`](gen_effect_calibration.md), where
individual_level effects of genotypes on taxa are computed.

## Usage

``` r
compute_beta_g(beta, genotypes, noise, taxa_scale)
```

## Arguments

- beta:

  A `matrix` of taxa (in rows) across SNPs (in columns) giving the
  multiplicative effect of genotype on taxa abundances. Typically
  obtained from
  [`compute_beta_matrix_cluster()`](compute_beta_matrix_cluster.md).

- genotypes:

  A `matrix` of genotypes with SNPs in rows and individuals in columns.

- noise:

  A `numeric` scalar indicating the standard deviation of the Gaussian
  noise to add.

- taxa_scale:

  A `numeric` scalar to scale the noise added to each OTUs abundance.

## Value

A `matrix` with taxa in roxs and individuals in columns, where each
element represents the abundance of a given taxa in an individual,
influenced by the genotype.

## See also

[`gen_effect_calibration()`](gen_effect_calibration.md)

## Examples

``` r
set.seed(123)
n_taxa <- 5
n_snp <- 10
n_ind <- 50
# Simulate a small beta matrix
beta <- matrix(sample(seq(0.1, 1, by = 0.05), size = n_taxa*n_snp, replace = TRUE), nrow = n_taxa, ncol = n_snp)

# Simulate a small genotype matrix
genotypes <- matrix(sample(0:2, size = n_snp*n_ind, replace = TRUE), nrow = n_snp, ncol = n_ind)

# Compute beta_g matrix with noise
beta_g <- RITHMS:::compute_beta_g(beta = beta,
                                  genotypes = genotypes,
                                  noise = 0.05,
                                  taxa_scale = 1)
```
