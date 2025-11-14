# Useful function to extract Shannon diversity from RITHMS output

The gets functions use the output of [`holo_simu()`](holo_simu.md) to
extract the information of interest from a given generation.
`get_mean_diversity`extract the average Shannon diversity from a
generation object.

## Usage

``` r
get_mean_diversity(data)
```

## Arguments

- data:

  List corresponding to one generation, as returned by
  [`holo_simu()`](holo_simu.md). Containing simulation output.

## Value

The average Shannon diversity within a given generation.

## Note

This function requires a previous call to
[`richness_from_abundances_gen()`](richness_from_abundances_gen.md) to
compute different types of diversity within the generation.

## See also

[`get_microbiomes()`](get_microbiomes.md),
[`richness_from_abundances_gen()`](richness_from_abundances_gen.md)

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
                              
# Extract Shannon diversity for each generations
## First step, compute richness from abundances
richness_from_abundances <- generations_simu[-1] %>% map(get_microbiomes) %>% map(richness_from_abundances_gen, size_rmultinom = 10000)
## size_rmultinom = 10000 according to DeruPops dataset

mean_shannon_diversity <- richness_from_abundances %>% map(get_mean_diversity)
} # }
if (FALSE) { # \dontrun{
library(magrittr)
library(purrr)

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

# To extract CLR abundances of the G1 generation microbiome
G1_microbiome <- get_microbiomes(generations_simu$G1, CLR = TRUE)

# To extract all the relative abundances of microbiomes from all generations
# substract metadata
microbiomes <- generations_simu[-1] %>% map(get_microbiomes)

# To extract mean phenotypes of all generations
mean_phenotypes <- generations_simu[-1] %>% map(get_mean_phenotypes)

# To extract phenotypic values for each individual of generation G1
G1_phenotypes_values <- get_phenotypes_value(generations_simu$G1)

# To extract the phenotype of each generation as the combined effects of microbiota and direct genetic effects
phenotypes <- generations_simu[-1] %>% map(get_phenotypes)

# To extract omega beta G values for each individual of generation G1
G1_om_beta_g <- get_om_beta_g(generations_simu$G1)

# To extract ID of selected individuals to each generation
selected_inds <- generations_simu[-1] %>% map(get_selected_ind) 

# To extract Shannon diversity
# First step, compute richness from abundances
richness_from_abundances <- generations_simu[-1] %>% map(get_microbiomes) %>% map(richness_from_abundances_gen)
shannon_diversity <- richness_from_abundances %>% map(get_mean_diversity)
} # }
```
