<!-- README.md is generated from README.Rmd. Please edit that file -->
inbreedR
========

![Build Status](https://travis-ci.org/mastoffel/inbreedR.svg?branch=master)

inbreedR provides functions and workflows for the analysis of inbreeding and heterozygosity-fitness correlations (HFCs) based on molecular markers such as microsatellites and SNPs. It has four main application areas:

-   Quantify variance in inbreeding level through estimation of identitiy disequilibria (g2), heterozygosity-heterozygosity correlations (HHC) and variance in standardized multilocus heterozygosity (sMLH)

-   Calculate g2 for large SNP datasets through the use of a formula that avoids double summations and has not been implemented in any software. The use of data.table and parallelization allow bootstrap and permutation tests with acceptable computation time

-   Estimate important parameters such as the influence of inbreeding on heterozygosity and fitness

-   Get insights on the sensitivity of these measures towards the number of genetic markers used in your study through re- and subsampling tests.

You can install:

-   the latest released version from CRAN (not yet) with

    ``` r
    install.packages("inbreedR")
    ```

-   the latest development version from github with

    ``` r
    if (packageVersion("devtools") < 1.6) {
      install.packages("devtools")
    }
    devtools::install_github("mastoffel/inbreedR", build_vignettes = TRUE)
    ```

If you encounter a clear bug or if you have any suggestions for improvement, just contact me: martin.adam.stoffel\[at\]gmail.com

Get started with inbreedR
-------------------------

To get started read the vignette:

``` r
vignette("inbreedR_step_by_step", package = "inbreedR")
```
