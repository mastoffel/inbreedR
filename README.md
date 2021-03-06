<!-- README.md is generated from README.Rmd. Please edit that file -->
inbreedR
========

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/inbreedR)](https://cran.r-project.org/package=inbreedR) ![Build Status](https://travis-ci.org/mastoffel/inbreedR.svg?branch=master) [![](http://cranlogs.r-pkg.org/badges/grand-total/inbreedR)](https://cran.r-project.org/package=inbreedR)

`inbreedR` provides functions and workflows for the analysis of inbreeding and heterozygosity-fitness correlations (HFCs) based on molecular markers such as microsatellites and SNPs. It has four main application areas:

-   Quantifying variance in inbreeding through estimation of identitiy disequilibria (g2), heterozygosity-heterozygosity correlations (HHC) and variance in standardized multilocus heterozygosity (sMLH)

-   Calculating g2 for small and large SNP datasets. The use of `data.table` and parallelization speed up bootstrapping and permutation tests

-   Estimating central parameters within HFC theory, such as the influence of inbreeding on heterozygosity and fitness, and their confidence intervals.

-   Exploring the sensitivity of these measures towards the number of genetic markers using simulations

You can install:

-   the latest released version from CRAN with

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

If you encounter bug or if you have any suggestions for improvement, just contact me: martin.adam.stoffel\[at\]gmail.com

Get started with inbreedR
-------------------------

To get started read the vignette:

``` r
vignette("inbreedR_step_by_step", package = "inbreedR")
```

### Citation

Stoffel, M. A., Esser, M., Kardos, M., Humble, E., Nichols, H., David, P., & Hoffman, J. I. (2016). inbreedR: an R package for the analysis of inbreeding based on genetic markers. Methods in Ecology and Evolution, 7(11), 1331-1339. <doi:10.1111/2041-210X.12588>
