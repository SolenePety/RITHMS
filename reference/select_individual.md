# Obtain id of the parents selected for the next generation based on the criteria chosen by the user

Obtain id of the parents selected for the next generation based on the
criteria chosen by the user

## Usage

``` r
select_individual(
  phenotypes,
  microbiomes,
  genotypes,
  beta,
  beta_otu,
  selection,
  size_selection_F,
  size_selection_M,
  selection_type,
  size_rmultinom,
  w.param
)
```

## Arguments

- phenotypes:

  Phenotype values of the current generation given as the result of the
  combined effects of the microbiota and direct genetic effects.
  Typically [`get_phenotypes()`](get_phenotypes.md) output or
  [`compute_phenotypes()`](compute_phenotypes.md).

- microbiomes:

  Abundances for each individual of the current generation. Note: It is
  necessary to transform the abundances if they are CLR-transformed with
  the following command. `microbiomes_clr |> t() |> clrInv() |> t()`

- genotypes:

  Genotypes values of the current generation. Typically from the MoBPS
  `get.geno()` function.

- beta:

  Beta matrix from
  [`compute_beta_matrix_cluster()`](compute_beta_matrix_cluster.md)
  output.

- beta_otu:

  Omega parameter from
  [`calibrate_params_phenotypes()`](calibrate_params_phenotypes.md)
  output.

- selection:

  bool, if selection process needed, DEFAULT = FALSE

- size_selection_F:

  percentage of female to select.

- size_selection_M:

  percentage of male to select.

- selection_type:

  mode of selection to be used, value in ("GB", "B", "G", "diversity",
  "div.GB"), DEFAULT = "GB"

- size_rmultinom:

  Integer; specifying the total number of object for the multinomial
  sampling(default: 10000, according to DeruPop.rds dataset).

- w.param:

  in case div.GB selection mode is chosen.
