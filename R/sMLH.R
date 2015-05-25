#' Calculate standardized multilocus heterozygosity (sMLH)
#' 
#' sMLH is defined as the total number of heterozygous loci in an individual divided 
#' by the sum of average observed heterozygosities in the population over the 
#' subset of loci successfully typed in the focal individual
#'
#' @param genotypes data.frame with individuals in rows and nloci in columns,
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
#' data(seal_microsats)
#' genotypes <- convert_raw(seal_microsats, miss = NA)
#' het <- sMLH(genotypes)
#'
#' @export
#'

sMLH <- function(genotypes) {
    genes <- as.matrix(genotypes)
    # get logical matrix of non-missing values as TRUE
    typed <- (genes == 1) | (genes == 0)
    # initialise
    nloc <- ncol(genes)
    nind <- nrow(genes)
    typed_sum <- colSums(typed)
    # heterozygosity per locus
    het_loc <- colSums(genes == 1) / typed_sum
    # replicate vector to matrix
    het_loc_mat <- matrix(het_loc, nrow = nind, ncol = nloc, byrow = TRUE)
    het_loc_mat[!typed] <- 0
    mh <- rowSums(het_loc_mat)
    N  <- rowSums(typed)
    H  <- rowSums(genes == 1)
    sMLH <- (H/N)/(mh/N)
    names(sMLH) <- row.names(genotypes)
    sMLH
}
