# g2_snps with data.table
library(data.table)
library(dplyr)
genotypes <- (tbl_df(fread("mice_snps.txt")))
genotypes[genotypes == -99] <- -1
genotypes[genotypes == NA] <- -1

genotypes <- read.table("raw_41loci_ordered.txt", row.names = 1)
genotypes <- convert_raw(genotypes, NAval = NA)

# genotypes2 <- tbl_df(read.table("mice_snps.txt"))

g2_snps_new <- function(genotypes) {
        # loci as rows
        origin <- t(genotypes)
        # define matrix with 1 for missing and 0 for all others
        m <- origin
        m[m==1] <- 0
        m[m==-1] <- 1
        # H matrix with 0 for -1
        h <- origin
        h[h==-1] <- 0
        # number of individuals/loci
        n <- ncol(origin)
        l <- nrow(origin)
        # vector with rowsums for missing data matrix
        m_loc_temp <- apply(m, 1, sum, na.rm = TRUE)
        m_loc <- m_loc_temp / n
        # rowsum h
        rowsum_h <- apply(h, 1, sum, na.rm = TRUE)
        colsum_h <- apply(h, 2, sum, na.rm = TRUE)
        # variance of colsums
        # B_hat <- var((apply(h, 2, sum, na.rm = TRUE)))
        # B_hat2 <- (1/(n-1)) * (sum(colsum_h^2) - (1/n) * sum(h)^2)
        # some sums
        h_sum <- sum(h, na.rm = TRUE)
        sum_rowsum_squared <- sum(rowsum_h^2, na.rm = TRUE)
        sum_colsum_squared <- sum(colsum_h^2, na.rm = TRUE)
        # C_hat <- (1/(n-1)) * (sum(h, na.rm = TRUE) - (1/n) * sum_rowsum_squared)
        # sum_colsum_squared <- sum(apply(h, 2, sum, na.rm = TRUE)^2)
        # A_hat <- (1/(n * (n-1))) * (sum(h, na.rm = TRUE)^2 - sum_rowsum_squared -
                                            # sum_colsum_squared - sum(h, na.rm = TRUE))

        # numerator
        numer <- (n-1) * (sum(colsum_h^2) - h_sum) /
                 (h_sum^2 - sum(rowsum_h^2) - sum(colsum_h^2) + h_sum)
        #
        M_ind <- vector()
        for (i in 1:n) {
                M_temp <- vector()
                for (k in 1:l) {
                        M_temp[k] <- (rowsum_h[k] * m[k, i]) / (1 - m_loc[k])
                }
                M_ind[i] <- (sum(M_temp, na.rm = TRUE))^2 - sum(M_temp^2, na.rm = TRUE)
        }
        M_hat <- (1/(n^2)) * sum(M_ind, na.rm = TRUE)

        #
        mat <- matrix(0, l, l)
        mat[upper.tri(mat)] <- 1

        #
        X_temp <- vector()
        for (k in 1:l) {
                X_temp[k] <- (rowsum_h[k] * m_loc[k]) / (1 - m_loc[k])
        }

        #
        a_hat_temp <- X_temp %*% (mat %*% matrix(X_temp, ncol = 1)) * (1/n) * 2

        # a_hat_temp <- (1/n) * (sum(X_temp^2, na.rm = TRUE) - (sum(X_temp, na.rm = TRUE))^2)

        a_hat_denom <- (1/n) * ((sum(h, na.rm = TRUE))^2 - sum_rowsum_squared)

        a_hat <- (M_hat - a_hat_temp) / a_hat_denom

        g2_snps <- (numer / (1 + a_hat)) - 1
}
