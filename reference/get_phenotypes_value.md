# Useful function to extract phenotype values from RITHMS output ("gq" by default)

The gets functions use the output of [`holo_simu()`](holo_simu.md) to
extract the information of interest from a given generation.
`get_phenotypes_values` extract the set of phenotype values from a
generation object.

## Usage

``` r
get_phenotypes_value(data, value = "gq")
```

## Arguments

- data:

  List corresponding to one generation, as returned by
  [`holo_simu()`](holo_simu.md). Containing simulation output.

- value:

  A String caracter that precise the object to extract from the
  generation object. Must be "gq" to extract phenotype values.

## Value

A `data.frame` of phenotype values for each individuals of a given
generation.

## See also

[`get_mean_phenotypes()`](get_mean_phenotypes.md),
[`get_microbiomes()`](get_microbiomes.md),
[`get_om_beta_g()`](get_om_beta_g.md),
[`get_selected_ind()`](get_selected_ind.md),
[`get_phenotypes()`](get_phenotypes.md)

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
                              
 #Extract phenotype values for each individuals of G1 generation
 G1_phenotypes_values <- get_phenotypes_value(generations_simu$G1)
 
 #Extract phenotype values for each individuals of all generations
 ## Don't forget to substract the metadata
 phenotypes_values <- generations_simu[-1] %>% map(get_phenotypes_value)
 } # }
 
```
