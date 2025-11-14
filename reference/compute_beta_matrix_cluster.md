# Generate beta matrix giving genetic effect per SNP on taxa abundances.

Generate beta matrix giving genetic effect per SNP on taxa abundances.

## Usage

``` r
compute_beta_matrix_cluster(
  n_b,
  n_g,
  n_clust,
  n_qtl_o,
  n_otus,
  effect_size = 1,
  correlation = 1,
  beta_info = NULL
)
```

## Arguments

- n_b:

  Number of taxa

- n_g:

  Number of SNPs

- n_clust:

  A vector with a length matching the total number of taxa with values
  from 0 to the number of clusters. Typically the
  [`assign_taxa()`](formatting_data.md) output.

- n_qtl_o:

  Number of causative QTL on taxa abundances (per taxon)

- n_otus:

  Number of taxa under genetic control.

- effect_size:

  Vector giving the size of genetic effect to try

- correlation:

  A numeric value between 0 and 1 representing the target correlation
  between taxa within the same cluster.

- beta_info:

  Informations from Beta matrix of genetic effects

## Value

A `list` of two objects linked to the beta matrix. The beta matrix
itself : a `matrix` of taxa (in rows) across SNPs (in columns) giving
the multiplicative effect of genotype on taxa abundances. A small
`data.frame` containing the parameters used during the beta matrix
simulation, including `cluster`, `id_otu` and `id_qtl_o`.

## See also

[`compute_beta_g()`](compute_beta_g.md), [`holo_simu()`](holo_simu.md)
