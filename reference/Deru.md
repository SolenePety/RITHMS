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

## Examples

``` r
# Access the dataset
data("Deru")

# View the first few rows of the taxa matrix
head(Deru$microbiome[1:5,1:5])
#>   OTU1 OTU2 OTU6793 OTU3 OTU4
#> 1   30  593       4  630  414
#> 2  254  275      62 1131  446
#> 3  181  487     164 1472  656
#> 4  469  665      45  640  328
#> 5  771  519      21  758  347

# Access the population object
population <- Deru$population
```
