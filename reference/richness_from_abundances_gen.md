# Estimate diversity metrics from relative abundances

This function estimate diversity metrics (Observed, Shannon, Inverse
Simpson) from the matrix of relative abundances (see
[`get_microbiomes()`](get_microbiomes.md)). It uses multinomial sampling
to simulate read counts from abundances, and computes diversity metrics
across `n_loop` in order to obtain robust estimation. This function is
particularly useful when selection is based on diversity.

## Usage

``` r
richness_from_abundances_gen(
  microbiome_matrix,
  size_rmultinom = 10000,
  n_loop = 10,
  plot = T
)
```

## Arguments

- microbiome_matrix:

  A matrix of relative abundances (individuals in rows and OTUs in
  columns, see [`get_microbiomes()`](get_microbiomes.md) output).

- size_rmultinom:

  Integer; specifying the total number of object for the multinomial
  sampling(default: 10000, according to DeruPop.rds dataset).

- n_loop:

  Integer; number of multinomial resampling iterations to perform
  (default: 10).

- plot:

  Logical; not currently implemented

## Value

A `data.frame`of average diversity metrics (Observed, Shannon, Inverse
Simpson) for each sample.

## See also

[`get_microbiomes()`](get_microbiomes.md),
[`phyloseq::estimate_richness()`](https://rdrr.io/pkg/phyloseq/man/estimate_richness.html)

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
                               
# Extract microbiomes matrix for each generations
microbiomes <- generations_simu[-1] %>% map(get_microbiomes)
 
# Estimate diversity metrics
richness_from_abundances <- microbiomes %>% map(richness_from_abundances_gen, size_rmultinom = 10000) 
## size_rmultinom = 10000 according to DeruPops dataset
} # }
```
