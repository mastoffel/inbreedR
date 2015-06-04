#' Genotype format converter
#'
#' Turns raw genotype data into 0 (homozygote), 1 (heterozygote) and NA (missing), which is the working format for 
#' the inbreedR functions.
#' A raw genotype matrix has individuals in rows and each locus in two adjacent columns. Individual ID´s can be rownames.
#' Type `data(mouse_msats)` for an example raw genotype data frame.
#'
#' @param genotypes Raw genotype data frame or matrix. Rows represent individuals and each locus has two adjacent columns. 
#'        Alleles within loci can be coded as numbers (e.g. microsatellite length) or characters (e.g. "A", "T")
#'        See data(mouse_msat) for an example. Missing values should be coded as NA or any negative number.
#' @param miss The value for missing data in the raw genotype data. Will be converted to NA.
#'
#' @return Data.frame object with 0 (homozygote), 1 (heterozygote) and NA (missing data).
#'         Each locus is thus represented by one column.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @examples
#' # Mouse microsatellite data with missing values coded as NA
#' data(mouse_msats)
#' genotypes <- convert_raw(mouse_msats, miss = NA)
#' head(genotypes)
#'
#' @export


convert_raw <- function(genotypes, miss = NULL) {
    
if (is.null(miss)) {
    stop("miss argument missing")
}

if (is.na(miss)) {
        genotypes[is.na(genotypes)] <- NA
} else if (is.numeric(miss)) {
        genotypes[genotypes == miss] <- NA
} else {
    stop("miss has to be NA or numeric")
}


check_het <- function(x) {

        if ((length(x) %% 2) == 1) {
                s1 <- seq(1, (length(x)-1), 2)
                newx <- as.vector(rep(NA, ceiling(length(x)/2)))
                }
        if ((length(x) %% 2) == 0) {
                s1 <- seq(1, length(x), 2)
                newx <- as.vector(rep(NA, length(x)/2))
                }

        count <- 1
        for(i in s1){
                if ((is.na(x[i]) | is.na(x[i + 1]))) {
                        newx[count] <- NA
                        count <- count + 1
                } else if (x[i] == x[i + 1]) {
                        newx[count] <- 0
                        count <- count + 1
                } else if (x[i] != x[i + 1]) {
                        newx[count] <- 1
                        count <- count + 1
                }
        }
        newx
}

# original full data matrix
origin <- as.data.frame(t(apply(genotypes, 1, check_het)))
row.names(origin) <- row.names(genotypes)
origin
}
