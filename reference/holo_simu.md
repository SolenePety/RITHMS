# Main function to generate transgenerational hologenomic data

Simulates hologenomic data over multiple generations with genetic,
microbial, and environmental effects.

## Usage

``` r
holo_simu(
  h2,
  b2,
  founder_object,
  n_ind = NULL,
  n_clust = NULL,
  n_gen = 5,
  qtn_y = NULL,
  correlation = 0.5,
  otu_g = 0.05,
  lambda = 0.5,
  effect.size = 0.1,
  mix.params = c(0.75, 0.25),
  mix.params.M = c(0.75, 0.25),
  noise.microbiome = 0.1,
  dir = T,
  ao = 25,
  size_rmultinom = 10000,
  selection = F,
  size_selection_F = NULL,
  size_selection_M = NULL,
  selection_type = "GB",
  w.param = c(0.5, 0.5),
  thetaX = NULL,
  env_gen = NULL,
  seed = 1234,
  verbose = T
)
```

## Arguments

- h2:

  direct heritability value, between 0 and 1.

- b2:

  microbiability value, between 0 and 1.

- founder_object:

  output of [`generate_founder()`](generate_founder.md) function.

- n_ind:

  number of individual per generation.

- n_clust:

  vector with taxa assignment, typical output of
  [`assign_taxa()`](formatting_data.md)

- n_gen:

  number of generation, DEFAULT = 5

- qtn_y:

  number of causal SNPs for the phenotypes.

- correlation:

  Correlation between taxa within the same cluster, value between 0 and
  1, DEFAULT = 0.5

- otu_g:

  percentage of taxa under genetic control, between 0 and 1, DEFAULT =
  0.05

- lambda:

  proportion of microbiome of offspring coming from vertical
  transmission, DEFAULT = 0.5

- effect.size:

  Vector giving the size of genetic effect to try.

- mix.params:

  Vector of two numeric values giving the weights for the regularisation
  of the base population microbiome. `mix.params[1]` = weight for raw
  microbiome and `mix.params[2]` = weight for mean microbiome. DEFAULT =
  c(0.75,0.25).

- mix.params.M:

  A vector of two numeric values specifying weights between Dirichlet
  samples and the original mean. DEFAULT = c(0.75,0.25). With
  `mix.params.M[1]` the dirichlet microbiome coefficient and
  `mix.params.M[2]` the mean microbiome coefficient.

- noise.microbiome:

  sd of microbiome noise, DEFAULT = 0.1

- dir:

  Logical; Mentions if the ambient microbiome is generated via a
  Dirichlet law or a only the `mean_microbiome`. DEFAULT = T

- ao:

  A numeric scalar used as the concentration parameter for the Dirichlet
  distribution.

- size_rmultinom:

  Integer; specifying the total number of object for the multinomial
  sampling(default: 10000, according to DeruPop.rds dataset).

- selection:

  bool, if selection process needed, DEFAULT = FALSE

- size_selection_F:

  percentage of female to select.

- size_selection_M:

  percentage of male to select.

- selection_type:

  mode of selection to be used, value in ("GB", "B", "G", "diversity",
  "div.GB"), DEFAULT = "GB"

- w.param:

  in case div.GB selection mode is chosen.

- thetaX:

  Optionnal matrix specifying environmental effects applied to the
  microbiome, such as antibiotic treatment. This matrix should be of
  dimension `n_taxa x n_individuals`, and typically constructed as the
  product of:

  - a vector of taxa-specific effects `theta` (negative or positive
    values for all taxa).

  - a binary vector `X` encoding individual exposure (1 for treated, 0
    ofr untreated). If `NULL` (default value), no environmental effect
    is applied. Example use case: refer to the vignette on [Generate
    figures](https://solenepety.github.io/RITHMS/articles/generate-figures.html#introduction-of-transient-perturbations-of-the-microbiota)

- env_gen:

  vector of booleans.

- seed:

  seed value for samplings in the function.

- verbose:

  bool, DEFAULT = T

## Value

An important list with several different objects :

- `parameters`: the set of parameters used for the simulation when
  calling `holo_simu()`

- `metadata`: metadata info such as beta matrix details

- `G0` to `G5` (5 generations by default): computational characteristics
  of generations. For each generation, genotypes, microbiomes,
  phenotypes, pedigree and individuals selected can be reachable.

## Required parameters

- `h2`

- `b2`

- `founder_object`

- `n_ind`

- `n_clust`

## Generation parameters

- `n_gen`

## Genotype-related parameters

- `qtn_y`

## Microbiome-related parameters

- `correlation`

- `otu_g`

- `lambda`

- `effect.size`

- `mix.params`

- `mix.params.M`

- `noise.microbiome`

- `dir`

- `ao`

## Selection-related parameters

- `selection`

- `size_selection_F`, `size_selection_M`

- `selection_type`

## Environmental effects

- `thetaX`

- `env_gen`

## Misc

- `seed`

- `verbose`

## See also

[`compute_beta_matrix_cluster()`](compute_beta_matrix_cluster.md),
[`compute_mean_microbiome()`](compute_mean_microbiome.md),
[`compute_current_microbiome()`](compute_current_microbiome.md),
[`compute_phenotypes()`](compute_phenotypes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  data("Deru")
  ToyData <- Deru
  taxa_assign_g <- assign_taxa(founder_object = ToyData)
  generations_simu <- holo_simu(h2 = 0.25,
                                 b2 = 0.25,
                                 founder_object = ToyData,
                                 n_clust = taxa_assign_g,
                                 n_ind = 500,
                                 verbose = FALSE,
                                 noise.microbiome = 0.5,
                                 effect.size = 0.3,
                                 lambda = 0.5,
                                 dir = TRUE,
                                 selection = FALSE,
                                 seed = 1234)
} # }  
```
