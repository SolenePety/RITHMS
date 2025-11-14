# Useful function to extract all phenotype values as data frame from RITHMS output

The gets functions use the output of [`holo_simu()`](holo_simu.md) to
extract the information of interest from a given generation.
`get_phenotypes`extract phenotype for each individuals as the combined
effects of microbiota and direct genetic effects, from a generation
object.

## Usage

``` r
get_phenotypes(data)
```

## Arguments

- data:

  List corresponding to one generation, as returned by
  [`holo_simu()`](holo_simu.md). Containing simulation output.

## Value

A `data.frame` of phenotype for each individuals of a given generation.
Phenotype values are given as the result of the combined effects of the
microbiota and direct genetic effects.

## See also

[`get_mean_phenotypes()`](get_mean_phenotypes.md),
[`get_phenotypes_value()`](get_phenotypes_value.md),
[`get_om_beta_g()`](get_om_beta_g.md),
[`get_selected_ind()`](get_selected_ind.md),
[`get_microbiomes()`](get_microbiomes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(magrittr)
library(purrr)
data("Deru")
ToyData <- Deru
taxa_assign_g <- assign_taxa(founder_object = ToyData)
generations_simu <- holo_simu(h2 = 0.25, b2 = 0.25, founder_object = ToyData,
                              n_clust = taxa_assign_g, n_ind = 500,
                              verbose = FALSE, seed = 1234)
                              
#Extract phenotypes of the G1 generation
G1_phenotypes <- get_phenotypes(generations_simu$G1)
 
#Extract phenotypes for all generations
## Don't forget to substract the metadata
phenotypes <- generations_simu[-1] %>% map(get_phenotypes)
} # }
```
