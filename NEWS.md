# inbreedR 0.2.0

## Improvements

* Bootstrapping over individuals for `r2_hf()` and `r2_Wf()`

* plotting histograms with CI for `r2_hf()` and `r2_Wf()`

* `r2_hf()` has an additional plot argument now, specify `plottype = "histogram"`to visualize
bootstrapping or `plottype = "boxplot"` to show the boxplots resulting from resampling of different
loci subsets.


# inbreedR 0.3.0

## Improvements

* `g2_resampling` function deleted 

* `simulate_g2` function added. This function simulates genotypes
from which different sized marker sets can be independently drawn to
estimate the precision and magnitude of g2 for a given dataset. Also works with larger
(SNP) datasets.

* `simulate_r2_hf` function added. This function uses the same simulation as `simulate_g2`
to estimate the expected correlation between heteorzygosity and inbreeding for 
varying number of markers. Also works with larger (SNP) datasets.

* `MLH` function added. MLH is the unstandardized version of the existing `sMLH` function.

* `subsets` argument in `r2_hf` function deprecated. Although you can infer the magnitude of
  the esimate by subsampling, the variation in estimates is biased. It is recommended to
  use the new `simulate_r2_hf` function instead.

# inbreedR 0.3.1

* deleted packages Hmisc and scales and exchanged with base R code

* authors added

* `r2_hf()` bound to 0-1

* verbose argument added to `g2_microsats` and `g2_snps`

* section on how to extract genotypes from VCF file added to vignette

# inbreedR 0.3.2

* added citation

* added section to vignette on calculating g2 in real-world SNP datasets

# inbreedR 0.3.3

* fixed linux CRAN built
