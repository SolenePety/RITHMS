
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RITHMS

<!-- badges: start -->
<!-- badges: end -->
<!-- Hexsticker -->

<img src="man/figures/hex_rithms.png" alt="Hexsticker" width="150" style="float: right; margin-left: 20px;" />

Our framework, an R Implementation of a Transgenerational Hologenomic
Model-based Simulator (RITHMS) an open-source package, builds upon the
MoBPS package and incorporates the distinctive architecture of the
microbiota, notably vertical and horizontal transfer as well as
modulation due to the environment and host genetics. In addition, RITHMS
can account for a variety of selection strategies, is adaptable to
different genetic architectures.

Full documentation website on: <https://SolenePety.github.io/RITHMS>

Here is a little summary of how RITHMS work, but you can read the
[preprint](https://github.com/SolenePety/RITHMS/blob/main/RITHMS_Pety.pdf)
for more details.

<figure>
<img src="man/figures/core_algorithm.png" alt="Schema" />
<figcaption aria-hidden="true">Schema</figcaption>
</figure>

## Installation

You can install the development version of RITHMS from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("SolenePety/RITHMS")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(RITHMS)
## basic example code
```

You already have a toy dataset that show you the expected structure of
input data :

``` r
# datafile <- system.file("BesePopTest.rds", package = "RITHMS")
# ToyData <- readRDS(datafile)
```

But we can do much more… you can check also [this
vignette](https://png.pngtree.com/png-vector/20220616/ourmid/pngtree-work-in-progress-warning-sign-with-yellow-and-black-stripes-painted-png-image_5060340.png)
to generate the beautiful figures coming from the article
