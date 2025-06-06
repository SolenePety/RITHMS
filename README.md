
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RITHMS <img src="man/figures/hex_rithms.png" align="right" width="155" height="180"/>

<!-- badges: start -->

Our framework, **R** **I**mplementation of a **T**ransgenerational
**H**ologenomic **M**odel-based **S**imulator (RITHMS) is an open-source
package, builds upon the MoBPS package and incorporates the distinctive
architecture of the microbiota, notably vertical and horizontal transfer
as well as modulation due to the environment and host genetics. In
addition, RITHMS can account for a variety of selection strategies, is
adaptable to different genetic architectures. <!-- badges: end -->

Full documentation website on: <https://SolenePety.github.io/RITHMS>  
Last code version on: <https://github.com/SolenePety/RITHMS>

Here is a little summary of how RITHMS work, but you can read the
[preprint](https://arxiv.org/abs/2502.07366) for more details.

![](man/figures/core_algorithm.png)

## Installation

You can install the development version of RITHMS from GitHub using
`remotes` as shown below.

``` r

# install.packages("remotes")

remotes::install_github("SolenePety/RITHMS")
```

## Toy dataset

You already have a toy dataset, a subset from [Déru et
al. 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7538339/), there are
**1845 species** and **780 individuals**, that show you the expected
structure of input data :

``` r
library(RITHMS)
data("Deru")
ToyData <- Deru
```

To import your own dataset, you can refer to the following
[vignette](https://solenepety.github.io/RITHMS/articles/import-data.html).

## Quick Start from the toy dataset

``` r
taxa_assign_g <- assign_taxa(ToyData)
generations_simu <- holo_simu(h2 = 0.25,
                              b2 = 0.25,
                              founder_object = ToyData,
                              n_clust = taxa_assign_g,
                              n_ind = 500)
```

## Quick Start from your own dataset

If you’re using your own dataset instead of the toy dataset, make sure
to calibrate the simulation parameters accordingly. In particular, the
`n_ind` (number of individuals in each generation), the genetic effect
size `effect.size` and the multinomial sampling size parameter
`size_rmultinom` used in holo_simu() should be consistent with the
characteristics of your dataset.

To help you with this, please refer to the [dedicated
vignette](https://solenepety.github.io/RITHMS/articles/calibrate-simulation-parameters.html).

``` r
# from VCf format
founder_object <- read_input_data(path_to_microbiome = "/path/to/microbiome.txt",
                                  path_to_pedmap = "/path/to/vcf/'prefix'",
                                  biome_id_column = "ind_id")
# from Ped/Map format
founder_object <- read_input_data(path_to_microbiome = "/path/to/microbiome.txt",
                                  path_to_pedmap = "/path/to/pedmap/'prefix'",
                                  biome_id_column = "ind_id")

taxa_assign_g <- assign_taxa(founder_object)
generations_simu <- holo_simu(h2 = 0.25,
                              b2 = 0.25, 
                              founder_object = founder_object,
                              n_clust = taxa_assign_g,
                              n_ind = 500,
                              effect.size = 0.3,
                              size_rmultinom = 10000)
# Choose n_ind, effect.size and size_rmultinom such that it is consistent with the initial dimensions of your data set
```

## To go further

If you’re interested into reproducing the figures of the article you can
take a look at [this
vignette](https://solenepety.github.io/RITHMS/articles/generate-figures.html#fine-selection-of-h2_d-b2-and-selection-schemes)
to generate the figures coming from the article
