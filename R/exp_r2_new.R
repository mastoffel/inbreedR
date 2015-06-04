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
#' data(mouse_msats)
#' genotypes <- convert_raw(mouse_msats, miss = NA)
#' (out <- exp_r2(genotypes, subsets = c(2,4,6,8,10,12), nboot = 1000, type = "msats"))
#' plot(out)
#' @export
#'
#'

exp_r2_new <- function(genotypes, subsets = NULL, nboot = 100, type = c("msats", "snps")) {
    
#     if (!(steps > 1) | (steps > ncol(genotypes))) {
#         stop("steps have to be at least two and smaller or equal than the number of markers used")
#     }
    gtypes <- as.matrix(genotypes)
    
    # check for subset sequence
    if ((sum(subsets < 2)) != 0) stop("You cannot subset less than 2 markers")
    if ((sum(subsets > ncol(genotypes)))!= 0) stop("The number of subsetted markers cannot exceed the overall number of markers")
    if (any(subsets%%1 != 0)) stop("All subsets have to be specified by integers")
    
    # sorting
    subsets <- sort(subsets)
    
    # check g2 function argument
    if (length(type) == 2){
        type <- "msats"
    } else if (!((type == "msats")|(type == "snps"))){
        stop("type argument needs to be msats or snps")
    } 
    
    # define calculation of expected r2
    if (type == "msats") {
        calc_r2 <- function(gtypes) {
            g2 <- g2_microsats(gtypes)[["g2"]]
            # according to the miller paper, negative g2´s are set to r2 = 0.
            if (g2 < 0) return(r2 <- 0)
            var_sh <- var(sMLH(gtypes))
            r2 <- g2/var_sh
            r2
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
    
#     # calculate sequence of loci numbers to draw
#     nloc_draw <- (round(ncol(genotypes)/ steps))
#     nloc_vec <- seq(from = nloc_draw, to = ncol(genotypes), by = nloc_draw)
    
    # initialise
    all_r2 <- matrix(nrow = nboot, ncol = length(subsets))
    
    calc_r2_sub <- function(gtypes, i) {
        loci <- sample((1:ncol(gtypes)), i)
        out <- calc_r2(gtypes[, loci])
        out
    }
    
    # counter
    step_num <- 1
    for (i in subsets) {
        
        cat("\n", "Iterating subset number ", step_num, " from ", length(subsets), sep = "")
        if (step_num == length(subsets)) {
            cat("\n", "Last subset!", sep = "")
        }
        
        all_r2[, step_num] <- replicate(nboot, calc_r2_sub(gtypes, i))
        
        step_num <- step_num + 1
    }
    
    # expected r2 per subset
    exp_r2_res <- data.frame(r2 = c(all_r2), nloc = factor(rep(subsets, each = nboot)))
    
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