# Useful function to extract omega beta G values from RITHMS output

The gets functions use the output of [`holo_simu()`](holo_simu.md) to
extract the information of interest from a given generation.
`get_om_beta_g`extract omega beta G values from a generation object.

## Usage

``` r
get_om_beta_g(data)
```

## Arguments

- data:

  List corresponding to one generation, as returned by
  [`holo_simu()`](holo_simu.md). Containing simulation output.

## Value

A `data.frame` of omega beta G values for each individuals of a given
generation.

## See also

[`get_mean_phenotypes()`](get_mean_phenotypes.md),
[`get_phenotypes_value()`](get_phenotypes_value.md),
[`get_microbiomes()`](get_microbiomes.md),
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
                              
#Extract omega beta G values for each individuals of G1 generation
G1_om_beta_g <- get_om_beta_g(generations_simu$G1)
 
#Extract omega beta G values for each individuals of all generations
## Don't forget to substract the metadata
om_beta_g <- generations_simu[-1] %>% map(get_om_beta_g)
} # }
 
```
