# Deru Dataset

A subset of data from [Déru et al.
2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7538339/) containing
hologenomic data from 750 pigs fed a conventional diet.

## Usage

``` r
Deru
```

## Format

A `list` containing two different matrices, one for the microbiome and
the other for the genotypes. `microbiome` : a matrix with 780 rows
(individuals) and 1845 taxa, where each element represents the abundance
of a specific taxon in a specific individual. `population` : object from
[MoBPS](https://pubmed.ncbi.nlm.nih.gov/32229505/) package which
provides the genotypes of each individual, encoded as 0, 1 ,2.

## Source

[Déru et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9194801/) for
microbiota data and [Déru et al.
2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7538339/) for genotypes
data.
