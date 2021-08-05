#' Helper function to estimating g2 from really large datasets
#'
#' @param genotypes \code{data.frame} with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote), 1 (heterozygote) and NA (missing)
#' @param nperm number or permutations for to estimate a p-value
#' @param nboot number of bootstraps to estimate a confidence interval
#' @param boot_over Bootstrap over individuals by specifying "inds" and over loci with "loci". Defaults to "ind".
#' @param CI confidence interval (default to 0.95)
#' @param parallel Default is FALSE. If TRUE, bootstrapping and permutation tests are parallelized 
#' @param ncores Specify number of cores to use for parallelization. By default,
#'        all available cores are used.
#' @param verbose If FALSE, nothing will be printed to show the status of bootstraps and permutations.
#' @keywords internal
#' 

subset_snps <- function(genotypes, nperm = 0, nboot = 0, boot_over = "inds", 
    CI = 0.95, parallel = FALSE, ncores = NULL, verbose = TRUE, subset_loci = 1000){
    
    
    g2_subset <- function(genotypes) {
        geno_sub <- genotypes[, sample(subset_loci, replace = FALSE)]
        g2_out <-  g2_snps(geno_sub, nperm = 0, nboot = 0, boot_over, CI, parallel, ncores, verbose)
    }
    
   
    
    
}