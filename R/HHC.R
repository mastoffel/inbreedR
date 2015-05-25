#' Calculates heterzygosity-heterozygosity correlations with 
#' standardized multilocus heterozygosities (sMLH)
#' 
#'
#' @param genotypes data.frame with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote) and 1 (heterozygote)
#' @param iter number of iterations, i.e. splittings of the dataset 
#'
#' @return
#' Vector of het-het correlations
#'
#' @references
#' Balloux, F., Amos, W., & Coulson, T. (2004). Does heterozygosity estimate inbreeding
#' in real populations?. Molecular Ecology, 13(10), 3021-3031.
#' 
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) 
#'        
#' @examples
#' data(seal_microsats)
#' genotypes <- convert_raw(seal_microsats, miss = NA)
#' out <- HHC(genotypes, iter = 100)
#'
#' @export
#'
#'

HHC <- function(genotypes, iter = 100) {
    # initialise
    loci <- 1:ncol(genotypes)
    loc_num <- ncol(genotypes)
    
    calc_cor <- function(num_iter, genotypes) {
        new_ord <- sample(loci)
        sMLH1 <- inbreedR::sMLH(genotypes[new_ord[1:floor(loc_num/2)]])
        sMLH2 <- inbreedR::sMLH(genotypes[new_ord[(floor(loc_num/2) + 1):loc_num]])
        het_het_cor <- cor(sMLH1, sMLH2)
        
        if (num_iter == 1) {
            cat("\n", "starting het-het correlations")
        } else if (num_iter == iter) {
            cat("\n", "done")
        } else if (num_iter %% 5 == 0) {
            cat("\n", num_iter, "iterations done")
        }
        return(het_het_cor)
    }
    het_het <- vapply(1:iter, calc_cor, numeric(1), genotypes)
#     out <- list(het_het = het_het, mean_het_het = mean(het_het))
#     class(het_het) <- "inbreedR"
    return(het_het)
}