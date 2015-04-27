#' Calculate standardized multilocus heterozygosity (sMLH) as well as a distribution
#' of sMLH under the hypothesis of no inbreeding (g2 = 0).
#' 
#' This function has the purpose of comparing the distribution of sMLH
#' in the sample to a theoretical distribution of sMLH when there is no inbreeding.
#' Distributions of sMLH under the hypothesis g2 = 0 are obtained by shuffling
#' genotypes at each locus randomly across individual.s
#'
#' @param genotypes data.frame with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote) and 1 (heterozygote)
#' @param nperm Number of permutations in order to get the distribution of sMLH values
#'        when there is no inbreeding.       
#'
#' @return
#' \item{sMLH_perm_mean}{Means of sMLH´s of permuted data.}
#' \item{sMLH_perm_mean}{Variances of sMLH´s of permuted data.}
#' \item{emp_sMLH}{empirical sMLH}
#' 
#' @references
#' Hoffman, J.I., Simpson, F., David, P., Rijks, J.M., Kuiken, T., Thorne, M.A.S., Lacey, R.C. & Dasmahapatra, K.K. (2014) High-throughput sequencing reveals inbreeding depression in a natural population.
#' Proceedings of the National Academy of Sciences of the United States of America, 111: 3775-3780. Doi: 10.1073/pnas.1318945111
#'
#'
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) 
#'        
#' @examples
#' data(seal_microsats)
#' # tranform raw genotypes into 0/1 format
#' genotypes <- convert_raw(seal_microsats, NAval = NA)
#' (het <- hetvar_nullhyp(genotypes, nperm = 100))
#'
#' @export

hetvar_nullhyp <- function(genotypes, nperm = 100) {
    # genotypes <- as.matrix(genotypes)
    
    permute <- function(genotypes) {
        geno_perm <- (apply(genotypes, 2, sample))
        perm_sMLH <- sMLH(geno_perm)
        sMLH_perm_mean <- mean(perm_sMLH)
        sMLH_perm_var <- var(perm_sMLH)
        out <- list(sMLH_perm_mean = sMLH_perm_mean, 
                    sMLH_perm_var  = sMLH_perm_var)
    }
    
    # permute genotypes
    all_permutes <- as.data.frame(replicate(nperm, permute(genotypes)))

    emp_sMLH <- sMLH(genotypes)
    
    res <- list(sMLH_perm_mean = unname(unlist(all_permutes[1, ])),
                sMLH_perm_var  = unname(unlist(all_permutes[2, ])),
                emp_sMLH = emp_sMLH)
    
    class(res) <- "inbreed"
    return(res)
}