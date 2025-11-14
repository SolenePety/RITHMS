# Compute phenotype values based on generated objects of the current generation

Compute phenotype values based on generated objects of the current
generation

## Usage

``` r
compute_phenotypes(
  X,
  B,
  otu_list,
  qtl_list,
  beta_qtl,
  beta_otu,
  Nqtl_y,
  Notu_y,
  se
)
```

## Arguments

- X:

  A `matrix` of the current generation genotypes, encoding with 0, 1 and
  2.

- B:

  A `matrix` of the CLR-transformed abundances of the current generation
  microbiome.

- otu_list:

  List of causal OTUs for the phenotypes.

- qtl_list:

  List of causal SNPs for the phenotypes.

- beta_qtl:

  Alpha, the regression coefficients corresponding to the QTL effects on
  the phenotype. See
  [`calibrate_params_phenotypes()`](calibrate_params_phenotypes.md).

- beta_otu:

  Omega, the regression coefficients corresponding to taxa effects on
  the phenotype. See
  [`calibrate_params_phenotypes()`](calibrate_params_phenotypes.md).

- Nqtl_y:

  Integer, number of causal SNPs for the phenotypes.

- Notu_y:

  Integer, number of causal OTUs for the phenotypes.

- se:

  Phenotypes Standard deviation, see
  [`calibrate_params_phenotypes()`](calibrate_params_phenotypes.md).

## See also

[`calibrate_params_phenotypes()`](calibrate_params_phenotypes.md)
