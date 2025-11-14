# Assign all taxa to a cluster, eventually under genetic control using hclust

Assign all taxa to a cluster, eventually under genetic control using
hclust

## Usage

``` r
assign_taxa(founder_object, taxa_g = 0.05)
```

## Arguments

- founder_object:

  Output of read_input_data function

- taxa_g:

  Percentage of taxa under genetic control, DEFAULT = 0.1

## Value

A vector with a length matching the total number of taxa with values
from 0 to the number of clusters, 0 corresponding to the non under
genetic control cluster

## Examples

``` r
data("Deru")
ToyData <- Deru
taxa_assign_g <- assign_taxa(founder_object = ToyData,
            taxa_g = 0.2)
```
