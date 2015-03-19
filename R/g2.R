#' estimating second-order heterozygosity from genotype data
#'
#' @param genotypes data.frame with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote) and 1 (heterozygote).
#' @param missing Value of missing data. -1 if data.frame is output from \link{convert_raw}.
#'
#' @return g2
#'
#' @references
#' DAVID, P., PUJOL, B., VIARD, F., CASTELLA, V. and GOUDET, J. (2007),
#' Reliable selfing rate estimates from imperfect population genetic data. Molecular Ecology,
#' 16: 2474–2487. doi: 10.1111/j.1365-294X.2007.03330.x
#'
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Mareike Esser (messer@@techfak.uni-bielefeld.de)
#'
#'
#' @examples
#' data(seal_microsats)
#' # tranform raw genotypes into 0/1 format
#' genotypes <- convert_raw(seal_microsats, NAval = NA)
#' g2_val <- g2(genotypes)
#'
#' @export
#'


g2 <- function(genotypes, missing = -1) {

        # transpose
        origin <- (t(genotypes))

        # check for any other missing or non 1, -1, 0 values and set them to missing (-1)
        all_out <- origin[!((origin == missing) | (origin == 0) | (origin == 1))]
        if (length(as.vector(all_out)) > 0) {
                warning(paste(length(all_out), "values found that do not equal 1, 0 or the defined missing value (", missing, ").
                              Those were set to missing."))
                origin[!((origin == missing) | (origin == 0) | (origin == 1))] <- -1
                origin[is.na(origin)] <- -1
        }

        # define matrix with 1 for missing and 0 for all others
        m <- origin
        m[m == 1] <- 0
        m[m == missing] <- 1

        # H matrix with 0 for missing
        h <- origin
        h[h == missing] <- 0

        n <- ncol(origin) # number of individuals
        l <- nrow(origin) # number of loci

        # mij: proportion of individuals missing on i and j ?s locus
        m_ij <- (m %*% t(m)) /n

        # vector with rowsums (proportion)
        m_loc <- rowSums(m) /n

        # numerator --------------------------------------------------------------------
        # pij entry says total amount of individuals that het locus i and locus j
        p <- h %*% t(h)
        missmat_num <- matrix(rep(0, l*l), ncol = l)

        # predefine vec
        vec <- c(1:nrow(h))


        for (i in seq(1:nrow(h))){
                for (j in seq(1:nrow(h))[-i]){
                        missmat_num[i,j]  <- 1/ (n * (1 - m_loc[i] - m_loc[j] + m_ij[i,j]))
                }
        }

#         for (i in seq(1:nrow(h))){
#                  vec <- vec[-i]
#                  missmat_num[vec, ] <- 1/ (n * (1 - m_loc[i] - m_loc[vec] + m_ij[vec, ]))
#         }

        numerator_mat <- p * missmat_num

       # rm(p)
        # rm(missmat_num)

        numerator <- sum(numerator_mat, na.rm = TRUE)

       # rm(numerator_mat)
        # denominator-------------------------------------------------------------------
        missmat_denom <- matrix(rep(0, l*l), ncol = l)

        # sum over loci

        for (i in seq(1:nrow(h))){
                for (j in seq(1:nrow(h))[-i]){
                        missmat_denom[i,j] <- 1/(((n * (n - 1) *
                                                      (1 - m_loc[i] - m_loc[j] + m_loc[i] * m_loc[j]))) -
                                                    (n * (m_ij[i, j] - m_loc[i] * m_loc[j])))

                }
        }

#         for (i in seq(1:nrow(h))){
#                 vec <- vec[-i]
#                 missmat_denom[i,vec] <- 1/((n - 1) * (1 - m_loc[i] - m_loc[vec] + m_loc[i] * m_loc[vec]) -
#                                                    (m_ij[i, vec] - m_loc[i] * m_loc[vec]))
#         }
        nullmat <- matrix(rep(1, n*n), ncol=n)
        diag(nullmat) <- 0
        q <- h %*% (nullmat %*% t(h))

        denominator_mat <-  missmat_denom * q
        denominator <- sum(denominator_mat, na.rm = TRUE)

        g2 <- (numerator / denominator) - 1

}
