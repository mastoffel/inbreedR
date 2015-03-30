#' estimating second-order heterozygosity from genotype data
#'
#' @param genotypes data.frame with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote) and 1 (heterozygote).
#' @param missing Value of missing data. -1 if data.frame is output from \link{convert_raw}.
#'
#' @return
#' \item{call}{Model call.}
#' \item{g2}{g2 value}
#' \item{p_val}{p value from permutation test}
#' \item{g2_permut}{g2 values from permuted genotypes}
#' \item{g2_boot}{g2 values from bootstrap samples}
#' \item{CI_boot}{confidence interval from bootstraps}
#' \item{se_boot}{standard error of g2 from bootstraps}
#' \item{nobs}{number of observations}
#' \item{nloc}{number of markers}
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
#' (g2_seals <- g2_microsats(genotypes, nperm = 100, nboot = 100, CI = 0.95))
#'
#' @export
#'


g2_microsats <- function(genotypes, nperm = 0, nboot = 0, CI = 0.95) {

        # transpose
        origin <- (t(genotypes))

        # check for any other missing or non 1, -1, 0 values and set them to missing (-1)
#         all_out <- origin[!((origin == missing) | (origin == 0) | (origin == 1))]
#         if (length(as.vector(all_out)) > 0) {
#                 warning(paste(length(all_out), "values found that do not equal 1, 0 or the defined missing value (", missing, ").
#                               Those were set to missing."))
#                 origin[!((origin == missing) | (origin == 0) | (origin == 1))] <- -1
#                 origin[is.na(origin)] <- -1
#         }

        # set all missings to -1
        origin[(origin!=0) & (origin!=1)] <- -1

        calc_g2 <- function(origin, perm = 1, boot = 1) {
                # define matrix with 1 for missing and 0 for all others
                m <- origin
                m[m == 1] <- 0
                m[m == -1] <- 1

                # H matrix with 0 for -1
                h <- origin
                h[(h!=0)&(h!=1)] <- 0

                n <- ncol(origin) # number of individuals
                l <- nrow(origin) # number of loci

                # mij: proportion of individuals missing on i and j ?s locus
                m_ij <- (m %*% t(m)) /n
                #diag(m_ij) <- 0
                # vector with rowsums (proportion)
                m_loc <- rowSums(m) /n

                # numerator --------------------------------------------------------------------
                # pij entry says total amount of individuals that are heterozygous at locus locus i and locus j
                p <- h %*% t(h)
                missmat_num <- matrix(rep(0, l*l), ncol = l)

                # predefine vec
                vec <- c(1:nrow(h))

                for (i in seq(1:nrow(h))){
                        vec_temp <- vec[-i]
                        missmat_num[i,  vec_temp] <- 1/ (n * (1 - m_loc[i] - m_loc[vec_temp] + m_ij[i,  vec_temp]))
                }

                numerator_mat <- p * missmat_num

                # rm(p)
                # rm(missmat_num)

                numerator <- sum(numerator_mat, na.rm = TRUE)

                # rm(numerator_mat)
                # denominator-------------------------------------------------------------------
                missmat_denom <- matrix(rep(0, l*l), ncol = l)

                for (i in seq(1:nrow(h))){
                        vec_temp <- vec[-i]
                        missmat_denom[i, vec_temp] <- 1/(n * (n - 1) * (1 - m_loc[i] - m_loc[vec_temp] + m_loc[i] * m_loc[vec_temp]) -
                                                                 (n * (m_ij[i, vec_temp] - m_loc[i] * m_loc[vec_temp])))
                }

                nullmat <- matrix(rep(1, n*n), ncol=n)
                diag(nullmat) <- 0
                q <- h %*% (nullmat %*% t(h))

                denominator_mat <-  missmat_denom * q
                denominator <- sum(denominator_mat, na.rm = TRUE)

                g2 <- (numerator / denominator) - 1

                if (perm %% 20 == 0) {
                        cat("\n", perm, "permutations done")
                } else if (perm == nperm-1) {
                        cat("\n", "### permutations finished ###")
                }

                if (boot %% 20 == 0) {
                        cat("\n", boot, "bootstraps done")
                } else if (boot == nboot-1) {
                        cat("\n", "### bootstrapping finished, hells yeah!! ###")
                }

                g2
        }

        # g2 point estimate
        g2_emp <- calc_g2(origin)

        # permutation of genotypes
        g2_permut <- rep(NA, nperm)
        p_permut <- NA

        if (nperm > 0) {
                #setkey(origin, eval(parse(names(origin)[1])))
                perm_genotypes <- function(perm, origin) {
                        # origin_perm <- origin[, lapply(.SD, sample)]
                        origin_perm <- t(apply(origin, 1, sample)) # to optimize
                        g2 <- calc_g2(origin_perm, perm = perm)

                }

                g2_permut <- c(g2_emp, sapply(1:(nperm-1), perm_genotypes, origin = origin))
                p_permut <- sum(g2_permut >= g2_emp) / nperm
                perm <- 1

        }
        g2_boot <- rep(NA, nboot)
        g2_se <- NA
        CI_boot <- c(NA,NA)

        if (nboot > 0) {

                boot_genotypes <- function(boot, origin) {
                        # bootstrap over individuals in columns
                        origin_boot <- origin[, sample(1:ncol(origin), replace = TRUE)]
                        g2 <- calc_g2(origin_boot, boot = boot)
                }

                g2_boot <- c(g2_emp, sapply(1:(nboot-1), boot_genotypes, origin = origin))
                g2_se <- sd(g2_boot)
                CI_boot <- quantile(g2_boot, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
        }

        res <- list(call=match.call(),
                    g2 = g2_emp, p_val = p_permut, g2_permut = g2_permut,
                    g2_boot = g2_boot, CI_boot = CI_boot, g2_se = g2_se,
                    nobs = nrow(genotypes), nloc = ncol(genotypes))

        class(res) <- "g2"
        return(res)

}
