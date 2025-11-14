# Useful function to extract ID of selected individuals from RITHMS output

The gets functions use the output of [`holo_simu()`](holo_simu.md) to
extract the information of interest from a given generation.
`get_selected_ind`extract selectd individuals IDs from a generation
object.

## Usage

``` r
get_selected_ind(data)
```

## Arguments

- data:

  List corresponding to one generation, as returned by
  [`holo_simu()`](holo_simu.md). Containing simulation output.

## Value

A`list` of the selected individuals of a given generation.

## See also

[`get_mean_phenotypes()`](get_mean_phenotypes.md),
[`get_phenotypes_value()`](get_phenotypes_value.md),
[`get_om_beta_g()`](get_om_beta_g.md),
[`get_microbiomes()`](get_microbiomes.md),
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
                              
#Extract selected individuals IDs for the G1 generation
G1_selected_ind <- get_selected_ind(generations_simu$G1)
 
#Extract selected individuals IDs for all generations
## Don't forget to substract the metadata
selected_ind <- generations_simu[-1] %>% map(get_selected_ind)
} # }
```
