#' Estimating g2 from SNP data
#'
#' @param genotypes data.frame with individuals in rows and loci in columns, containing genotypes coded as 0 (homozygote) or 1 (heterozygote)
#' @param nperm number or permutations for calculating the p-value
#' @param nboot number of bootstraps for CI
#' @param CI confidence interval
#'
#' @details Calculates g2 from big SNP datasets. Use convert_raw to convert raw genotypes (with 2 columns per locus) into
#'          the required format
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
#' Hoffman, J.I., Simpson, F., David, P., Rijks, J.M., Kuiken, T., Thorne, M.A.S., Lacey, R.C. & Dasmahapatra, K.K. (2014) High-throughput sequencing reveals inbreeding depression in a natural population.
#' Proceedings of the National Academy of Sciences of the United States of America, 111: 3775-3780. Doi: 10.1073/pnas.1318945111
#'
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Mareike Esser (messer@@techfak.uni-bielefeld.de)
#'
#' @examples
#' # load SNP genotypes in 0 (homozygous), 1(heterozygous) format.
#'
#' data(mice_snp_genotypes)
#' (g2_val <- g2_snps(mice_snp_genotypes, nperm = 100, nboot = 100, CI = 0.95))
#'
#'
#' @export

g2_snps <- function(genotypes, nperm = 0, nboot = 0, CI = 0.95) { # , missing = -1
        # transpose for congruency with formulae in paper
        origin <- t(genotypes)
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

                # vector with rowsums for missing data matrix
                m_loc_temp <- rowSums(m, na.rm = TRUE)
                m_loc <- m_loc_temp / n

                # precalculating sums
                rowsum_h <- rowSums(h, na.rm = TRUE)
                colsum_h <- colSums(h, na.rm = TRUE)
                sum_rowsum_squared <- sum(rowsum_h^2)
                sum_colsum_squared <- sum(colsum_h^2)
                h_sum <- sum(colsum_h, na.rm = TRUE)

                # numerator
                numer <- (n-1) * (sum_colsum_squared - h_sum) /
                        (h_sum^2 - sum_rowsum_squared - sum_colsum_squared + h_sum)

                M_temp <- rowsum_h * m / (1 - m_loc)
                M_ind <- colSums(M_temp, na.rm = TRUE)^2 - colSums(M_temp^2, na.rm = TRUE)

                M_hat <- (1/(n)) * sum(M_ind, na.rm = TRUE)
                X_temp <- rowsum_h * m_loc / (1-m_loc)

                a_hat_temp <- (sum(X_temp^2, na.rm = TRUE) - (sum(X_temp, na.rm = TRUE))^2)
                a_hat <- (M_hat + a_hat_temp) / (h_sum^2 - sum_rowsum_squared)
                g2_emp <- numer / (1 + a_hat) - 1

                if (perm %% 5 == 0) {
                        cat("\n", perm, "permutations done")
                } else if (perm == nperm-1) {
                        cat("\n", "### permutations finished ###")
                }

                if (boot %% 5 == 0) {
                        cat("\n", boot, "bootstraps done")
                } else if (boot == nboot-1) {
                        cat("\n", "### bootstrapping finished ###")
                }

                g2_emp
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
                        g2 <- calc_g2(origin_perm, perm)
                }

                g2_permut <- c(g2_emp, sapply(1:(nperm-1), perm_genotypes, origin))
                p_permut <- sum(c(g2_emp, g2_permut) >= g2_emp) / nperm
                perm <- 1

        }

        g2_boot <- rep(NA, nboot)
        g2_se <- NA
        CI_boot <- c(NA,NA)

        if (nboot > 0) {

                boot_genotypes <- function(boot, origin) {
                        # bootstrap over individuals in columns
                        origin_boot <- origin[, sample(1:ncol(origin), replace = TRUE)]
                        g2 <- calc_g2(origin_boot, perm, boot)
                }

                g2_boot <- c(g2_emp, sapply(1:(nboot-1), boot_genotypes, origin))
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
