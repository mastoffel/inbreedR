<!-- README.md is generated from README.Rmd. Please edit that file -->

# inbreedR

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/inbreedR)](https://cran.r-project.org/package=inbreedR)
![Build
Status](https://travis-ci.org/mastoffel/inbreedR.svg?branch=master)
[![](http://cranlogs.r-pkg.org/badges/grand-total/inbreedR)](https://cran.r-project.org/package=inbreedR)
[![Codecov test
coverage](https://codecov.io/gh/mastoffel/inbreedR/branch/master/graph/badge.svg)](https://codecov.io/gh/mastoffel/inbreedR?branch=master)

### Goal

`inbreedR` provides functions and workflows for the analysis of
inbreeding and heterozygosity-fitness correlations (HFCs) based on
molecular markers such as microsatellites and SNPs. In case of genomic
data, it’s most useful for lower density datasets where it is unclear
whether genotyped markers represent genome-wide diversity / inbreeding.
It has four main application areas:

-   Quantifying variance in inbreeding through estimation of identitiy
    disequilibria (g2), heterozygosity-heterozygosity correlations (HHC)
    and variance in standardized multilocus heterozygosity (sMLH)

-   Calculating g2 for small and large SNP datasets. The use of
    `data.table` and parallelization speed up bootstrapping and
    permutation tests

-   Estimating central parameters within HFC theory, such as the
    influence of inbreeding on heterozygosity and fitness, and their
    confidence intervals.

-   Exploring the sensitivity of these measures towards the number of
    genetic markers using simulations

### Installation

You can install the stable version of `inbreedR` from CRAN with:

``` r
install.packages("rptR")
```

Or the development version from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("mastoffel/inbreedR", build_vignettes = TRUE, dependencies = TRUE) 
# manual
browseVignettes("inbreedR")
```

If you find a bug, please report a minimal reproducible example in the
[issues](https://github.com/mastoffel/inbreedR/issues).

## Get started with inbreedR

To get started read the vignette:

``` r
vignette("inbreedR_step_by_step", package = "inbreedR")
```

### Citation

Stoffel, M. A., Esser, M., Kardos, M., Humble, E., Nichols, H., David,
P., & Hoffman, J. I. (2016). inbreedR: an R package for the analysis of
inbreeding based on genetic markers. *Methods in Ecology and Evolution*,
**7**(11), 1331-1339. <doi:10.1111/2041-210X.12588>
