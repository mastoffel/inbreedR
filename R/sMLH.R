
#' Calculate standardized multilocus heterozygosity (sMLH)
#' 
#' sMLH is defined as total number of heterozygous loci in an individual divided 
#' by the sum of average observed heterozygosities in the population over the 
#' subset of loci successfully typed in the focal individual
#'
#' @param genotypes data.frame with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote) and 1 (heterozygote)
#'
#' @return
#' Vector of individual standardized multilocus heterozygosities
#'
#' @references
#' Coltman, D. W., Pilkington, J. G., Smith, J. A., & Pemberton, J. M. (1999). 
#' Parasite-mediated selection against inbred Soay sheep in a free-living, 
#' island population. Evolution, 1259-1267.
#'
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) 
#'        
#' @examples
#' data(mice_snp_genotypes)
#' het <- sMLH(mice_snp_genotypes)
#'
#' @export
#'

sMLH <- function(genotypes) {
    
    genotypes <- as.matrix(genotypes)
    # get logical matrix of non-missing values as TRUE
    typed <- (genotypes == 1) | (genotypes == 0)
    # initialise
    loc <- ncol(genotypes)
    indiv <- nrow(genotypes)
    het_loc <- rep(NA, loc)
    typed_sum <- colSums(typed)
    # at each locus, which percent of individuals is het (without missings)
    het_loc <- colSums(genotypes == 1) / typed_sum
    het_loc_mat <- matrix(het_loc, nrow = indiv, ncol = loc, byrow = TRUE)
    
    het_loc_mat[!typed] <- 0
    mh <- rowSums(het_loc_mat)
    N  <- rowSums(typed)
    H  <- rowSums(genotypes == 1)
    sMLH <- (H/N)/(mh/N)
    
    names(sMLH) <- row.names(genotypes)
    sMLH
}
