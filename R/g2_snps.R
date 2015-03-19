#' Estimating g2 from SNP data
#'
#' @param genotypes data.frame with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote) or 1 (heterozygote)
#'
#' @details Calculates g2 from big SNP datasets. Use convert_raw to convert raw genotypes (with 2 columns per locus) into
#'          the required format
#'
#' @return g2 value
#'
#' @references
#' Hoffman, J.I., Simpson, F., David, P., Rijks, J.M., Kuiken, T., Thorne, M.A.S., Lacey, R.C. & Dasmahapatra, K.K. (2014) High-throughput sequencing reveals inbreeding depression in a natural population.
#' Proceedings of the National Academy of Sciences of the United States of America, 111: 3775-3780. Doi: 10.1073/pnas.1318945111
#'
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Mareike Esser (messer@@techfak.uni-bielefeld.de)
#'
#' @examples
#' # load SNP genotypes, already in 0 (homozygous), 1(heterozygous) format.
#'
#' data(mice_snp_genotypes)
#' g2_val <- g2_snps(genotypes)
#'
#'
#' @export

g2_snps <- function(genotypes) {

        origin <- t(genotypes)
        # define matrix with 1 for missing and 0 for all others
        m <- origin
        m[m==1] <- 0
        m[m==-1] <- 1

        # H matrix with 0 for -1
        h <- origin
        h[h==-1] <- 0

        n <- ncol(origin) # number of individuals
        l <- nrow(origin) # number of loci

        # vector with rowsums for missing data matrix
        m_loc_temp <- rowSums(m, na.rm = TRUE)
        m_loc <- m_loc_temp / n

        # precalculating sums
        rowsum_h <- rowSums(h, na.rm = TRUE)
        colsum_h <- colSums(h, na.rm = TRUE)
        sum_rowsum_squared <- sum(rowsum_h^2)
        sum_colsum_squared <- sum(colsum_h^2)
        h_sum <- sum(h, na.rm = TRUE)

        # numerator
        numer <- (n-1) * (sum_colsum_squared - h_sum) /
                 (h_sum^2 - sum_rowsum_squared - sum_colsum_squared + h_sum)

        M_temp <- rowsum_h * m / (1 - m_loc)
        M_ind <- colSums(M_temp, na.rm = TRUE)^2 - colSums(M_temp^2, na.rm = TRUE)

        M_hat <- (1/(n)) * sum(M_ind, na.rm = TRUE)
        X_temp <- rowsum_h * m_loc / (1-m_loc)

        a_hat_temp <- (sum(X_temp^2, na.rm = TRUE) - (sum(X_temp, na.rm = TRUE))^2)
        a_hat <- (M_hat + a_hat_temp) / (h_sum^2 - sum_rowsum_squared)
        g2_val <- numer / (1 + a_hat) - 1

}
