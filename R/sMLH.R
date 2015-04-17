# total number
# of heterozygous loci in an individual divided by the sum of average observed
# heterozygosities in the population over the subset of loci successfully typed in
# the focal individual

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
#' data(seal_microsats)
#' # tranform raw genotypes into 0/1 format
#' genotypes <- convert_raw(seal_microsats, NAval = NA)
#' het <- sMLH(genotypes)
#'
#' @export
#'

sMLH <- function(genotypes) {
    
    genotypes <- as.matrix(genotypes)
    typed <- (genotypes == 1) | (genotypes == 0)
    loc <- ncol(genotypes)
    indiv <- nrow(genotypes)
    het_loc <- vector()
    typed_sum <- colSums(typed)
    # at each locus, which percent of individuals is het (without missings)
    
    for (l in 1:loc) {
        het_loc[l] <- sum(genotypes[, l] == 1) / typed_sum[l]
    }
    
    sMLH <- rep(NA, nrow(genotypes))
    for (i in 1:indiv) {
        mh <- sum(het_loc[typed[i, ]])
        N <- sum(typed[i, ])
        H <- sum(genotypes[i, ] == 1)
        sMLH[i] <- (H/N)/(mh/N)
    }
    names(sMLH) <- rownames(genotypes)
    sMLH
}



# sMLH <- function (genotypes) {
#     genotypes <- as.matrix(genotypes)
#     typed <- (genotypes == 1) | (genotypes == 0)
#     
#     indiv <- nrow(genotypes)
#     all_typed_het <- vector()
#     
#     for (i in 1:indiv) {
#         all_typed_het[i] <- mean(rowSums(genotypes[, typed[i, ]]))
#         # apply(genotypes[, typed[i, ]], 1, function(x) sum(x[x==1])))
#     }
#     
#     het_indiv <- apply(genotypes, 1, function(x) sum(x[x==1]))
#     out <- het_indiv / all_typed_het
#     # logical matrix
#     out
# }