#' Expected r2 between sMLH and inbreeding level
#' 
#' Currently just working for microsat data!
#'
#' @param genotypes data.frame with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote) and 1 (heterozygote)
#' @param steps specifies number of subsets for increasing marker number
#' @param nboot number of bootstraps per step
#' @param type specifies g2 formula to take. Type "snps" for large datasets and "msats" for smaller datasets.
#' 
#' @return 
#' \item{call}{function call.}
#' \item{exp_r2_res}{expected r2 for each randomly subsetted dataset}
#' \item{summary_exp_r2}{r2 mean and sd for each number of subsetted loci}
#' \item{nobs}{number of observations}
#' \item{nloc}{number of markers}
#' 
#' @references
#' Szulkin, M., Bierne, N., & David, P. (2010). HETEROZYGOSITY-FITNESS CORRELATIONS: A TIME FOR REAPPRAISAL. 
#' Evolution, 64(5), 1202-1217.
#' 
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) 
#'        
#' @examples
#' data(seal_snps)
#' (out <- exp_r2(seal_snps, steps = 5, nboot = 5, type = "snps"))
#' plot(out)
#' @export
#'
#'

exp_r2 <- function(genotypes, steps = 10, nboot = 100, type = c("msats", "snps")) {
    
    if (!(steps > 1) | (steps > ncol(genotypes))) {
        stop("steps have to be at least two and smaller or equal than the number of markers used")
    }
    
    gtypes <- as.matrix(genotypes)
    
    if (length(type) == 2){
        type <- "msats"
    } else if (!((type == "msats")|(type == "snps"))){
        stop("type argument needs to be msats or snps")
    } 
    
    if (type == "msats") {
        calc_r2 <- function(gtypes) {
            g2 <- g2_microsats(gtypes)[["g2"]]
            # according to the miller paper, negative g2´s are set to r2 = 0.
            if (g2 < 0) return(r2 <- 0)
            var_sh <- var(sMLH(gtypes))
            r2 <- g2/var_sh
        }
    }
    
    if (type == "snps") {
        calc_r2 <- function(gtypes) {
            g2 <- g2_snps(gtypes)[["g2"]]
            # according to the miller paper, negative g2´s are set to r2 = 0.
            if (g2 < 0) return(r2 <- 0)
            var_sh <- var(sMLH(gtypes))
            r2 <- g2/var_sh
        }
    }
    
    # calculate sequence of loci numbers to draw
    nloc_draw <- (round(ncol(genotypes)/ steps))
    nloc_vec <- seq(from = nloc_draw, to = ncol(genotypes), by = nloc_draw)
    # initialise
    all_r2 <- matrix(nrow = nboot, ncol = length(nloc_vec))
    
    calc_r2_sub <- function(gtypes, i) {
        loci <- sample((1:ncol(gtypes)), i)
        out <- calc_r2(gtypes[, loci])
    }
    
    # counter
    step_num <- 1
    for (i in seq(1:length(nloc_vec))) {
        
        cat("\n", "Iterating subset number ", step_num, " from ", steps, sep = "")
        if (step_num == steps) {
        cat("\n", "Last subset!", sep = "")
        }
        
        all_r2[, i] <- replicate(nboot, calc_r2_sub(gtypes, nloc_vec[i]))
        
        step_num <- step_num + 1
    }
    
    # expected r2 per subset
    exp_r2_res <- data.frame(r2 = c(all_r2), nloc = factor(rep(nloc_vec, each = nboot)))
    
    # mean and sd per number of loci
    summary_exp_r2 <- aggregate(r2 ~ nloc, data = exp_r2_res, 
                             FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)))
    
    res <- list(call = match.call(),
                exp_r2_res = exp_r2_res,
                summary_exp_r2 = summary_exp_r2,
                nobs = nrow(genotypes), 
                nloc = ncol(genotypes))
    
    class(res) <- "inbreed"
    return(res)
    
}