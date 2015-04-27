sMLH2 <- function(genotypes) {
    
    genotypes <- as.data.frame(genotypes)
    typed <- (genotypes == 1) | (genotypes == 0)
    loc <- ncol(genotypes)
    indiv <- nrow(genotypes)
    het_loc <- vector()
    typed_sum <- colSums(typed)
    # at each locus, which percent of individuals is het (without missings)
    
         for (l in 1:loc) {
             het_loc[l] <- sum(genotypes[, l] == 1) / typed_sum[l]
         }
    #     
    # het_loc <- colSums(genotypes == 1) / typed_sum
    
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