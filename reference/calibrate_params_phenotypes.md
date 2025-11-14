# Calibration coefficients to compute phenotypes

This function is part of the first step of
[`holo_simu()`](holo_simu.md), focused on calibrating phynotypes
parameters before looping on the generations. All coefficents are
rescale during the process to ensure that standard deviation = 1 and the
variance of microbiote effect and genetic effect satisfy the target
microbiability and direct heritability.

## Usage

``` r
calibrate_params_phenotypes(
  X0,
  B0,
  h2,
  b2,
  otu_list,
  Nqtl_y,
  Notu_y = length(otu_list)
)
```

## Arguments

- X0:

  Matrix of the primary generation genotypes, given from
  [get.geno()](https://rdrr.io/pkg/MoBPS/man/get.geno.html) MoBPS
  function. SNPs are in rows and individuals in columns.

- B0:

  Matrix of the primary generation microbiomes, as CLR abundances. OTUs
  are in rows and individuals in columns.

- h2:

  direct heritability value, between 0 and 1.

- b2:

  microbiability value, between 0 and 1.

- otu_list:

  List of causal OTUs for the phenotypes.

- Nqtl_y:

  Integer; number of causal SNPs for the phenotypes.

- Notu_y:

  Integer; number of causal OTUs for the phenotypes.

## Value

A `list` of phenotype parameters such as alpha, omega, list of causal
SNPs for the phenotypes, standard deviation (=1) of phenotypes

## See also

[`compute_phenotypes()`](compute_phenotypes.md)

## Examples

``` r
library(magrittr)
n_ind <- 10
n_snp <- 50
n_otu <- 50

# Simulate a small genotype matrix 
set.seed(123)
X0 <- matrix(sample(0:2, n_snp * n_ind, replace = TRUE), nrow = n_snp, ncol = n_ind)

# Simulate a small microbiome counts matrix
B0_counts_table <- matrix(abs(rnorm(n_snp * n_ind, mean = 10, sd = 3)), nrow = n_otu, ncol = n_ind)

# Transform to relative abundances per individuals and then apply CLR transformation
B0_abund <- apply(B0_counts_table, 2, function(x) x/sum(x))
B0 <- compositions::clr(t(B0_abund)) %>% t()

# Randomly select causal OTUs
rownames(B0) <- paste0("OTU", 1:n_otu)
otu_list <- c(y = sample(rownames(B0), 10))


params <- suppressWarnings(RITHMS:::calibrate_params_phenotypes(X0 = X0,
                                      B0 = B0,
                                      h2 = 0.25,
                                      b2 = 0.25,
                                      otu_list = otu_list,
                                      Nqtl_y = 10,
                                      Notu_y = length(otu_list)))
str(params)
#> List of 4
#>  $ alpha   : num [1:10] 0.000119 -0.139382 0.718991 -0.084337 -0.217468 ...
#>  $ omega   : num [1:10] 1.147 1.564 1.253 -0.885 0.279 ...
#>  $ qtl_list: int [1:10] 28 15 11 29 17 41 25 12 19 7
#>  $ se      : num [1, 1] 1
 
```
