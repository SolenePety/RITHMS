# Useful function to extract mean phenotype values from RITHMS output

The gets functions use the output of [`holo_simu()`](holo_simu.md) to
extract the information of interest from a given generation.
`get_mean_phenotypes`extract average phenotype from a generation object.

## Usage

``` r
get_mean_phenotypes(data)
```

## Arguments

- data:

  List corresponding to one generation, as returned by
  [`holo_simu()`](holo_simu.md). Containing simulation output.

## Value

The mean phenotype value for a given generation.

## See also

[`get_microbiomes()`](get_microbiomes.md),
[`get_phenotypes_value()`](get_phenotypes_value.md),
[`get_om_beta_g()`](get_om_beta_g.md),
[`get_selected_ind()`](get_selected_ind.md),
[`get_phenotypes()`](get_phenotypes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(purrr)
library(magrittr)
data("Deru")
ToyData <- Deru
taxa_assign_g <- assign_taxa(founder_object = ToyData)
generations_simu <- holo_simu(h2 = 0.25, b2 = 0.25, founder_object = ToyData,
                              n_clust = taxa_assign_g, n_ind = 500,
                              verbose = FALSE, seed = 1234)
                              
# Extract mean phenotype value for G1 generation
G1_mean_phenotype <- get_mean_phenotypes(generations_simu$G1)

# Extract mean phenotype values of all generations
## Don't forget to substract the metadata
mean_phenotypes <- generations_simu[-1] %>% map(get_mean_phenotypes)
} # }
```
