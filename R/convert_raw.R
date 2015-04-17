#' Genotype format converter
#'
#' Turns raw genotype data into 0 (homozygote), 1 (heterozygote) and -1 (NA on at least one locus).
#' A raw genotype matrix has individuals in rows and loci in columns, while each locus has two adjacent.
#' Type \emph{data(seal_microsats)} for an example data frame.
#'
#' @param genotypes Raw genotype data frame or matrix. Has individuals in rows and each locus in two columns
#' @param NAval The prespecified value for missing elements in the data frame. Will be converted to -1 in the output.
#'
#'
#' @return Data.frame object with 0 (homozygote), 1 (heterozygote) and -1 (NA on at least one locus).
#'         Each locus is thus represented by one column.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @examples
#' # Seal microsatellite data with missing values specified with -99
#' data(seal_microsats)
#' genotypes <- convert_raw(seal_microsats, NAval = NA)
#' head(genotypes)
#'
#' @export



convert_raw <- function(genotypes, NAval = -99) {

if (is.na(NAval)) {
        genotypes[is.na(genotypes)] <- NA
} else {
        genotypes[genotypes == NAval] <- NA
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
                if (is.na(x[i] | x[i + 1])) {
                        newx[count] <- -1
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
