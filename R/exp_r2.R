#' Expected r2 between sMLH and inbreeding level
#' 
#'
#' @param genotypes data.frame with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote) and 1 (heterozygote)
#'
#'
#' @references
#' Szulkin, M., Bierne, N., & David, P. (2010). HETEROZYGOSITY-FITNESS CORRELATIONS: A TIME FOR REAPPRAISAL. 
#' Evolution, 64(5), 1202-1217.
#' 
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) 
#'        
#' @examples
#' data(mice_snp_genotypes)
#' result <- exp_r2(mice_snp_genotypes)
#' 
#' data(seal_microsats)
#' genotypes <- convert_raw(seal_microsats, miss_val = NA)
#'
#' @export
#'
#'

exp_r2 <- function(genotypes, parts = 10, nboot = 1000) {
    
    gtypes <- as.matrix(genotypes)
    
    calc_r2 <- function(gtypes) {
    g2 <- g2_microsats(gtypes)[["g2"]]
    var_sh <- var(sMLH(gtypes))
    r2 <- g2/var_sh
    }
    
    # calculate sequence of loci numbers to draw
    nloc_draw <- (floor(ncol(genotypes)/ parts))
    nloc_vec <- seq(from = nloc_draw, to = ncol(genotypes), by = nloc_draw)
    
    all_r2 <- matrix(nrow = nboot, ncol = parts)
    
    calc_r2_sub <- function(gtypes, i) {
        loci <- sample((1:ncol(gtypes)), i)
        out <- calc_r2(gtypes[, loci])
    }
    
    for (i in seq(1:parts)) {
        all_r2[, i] <- replicate(nboot, calc_r2_sub(gtypes, nloc_vec[i]))
    }
    
    res <- data.frame(r2 = c(all_r2), npart = factor(rep(1:parts, each = nboot)))
    
}