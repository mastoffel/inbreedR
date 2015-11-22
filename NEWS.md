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



