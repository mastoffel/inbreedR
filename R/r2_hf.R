#' Expected r2 between standardized multilocus heterozygosity (h) and inbreeding level (f)
#' 
#' 
#'
#' @param genotypes data.frame with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote) and 1 (heterozygote)
#' @param subsets a vector specifying the sizes of subsets to draw. For a subset of 20 markers, subsets = c(2, 5, 10, 15, 20) could
#'        be a reasonable choice. The minimum subset size is 2 and the maximum is the number of markers in the data.
#' @param nboot number re-draws per subset size.
#' @param type specifies g2 formula to take. Type "snps" for large datasets and "msats" for smaller datasets.
#' @param parallel Default is FALSE. If TRUE, bootstrapping and permutation tests are parallelized 
#' @param ncores Specify number of cores to use for parallelization. By default,
#'        all available cores are used.
#'        
#' 
#' @return 
#' \item{call}{function call.}
#' \item{r2_hf_full}{expected r2 between inbreeding and sMLH for the full dataset}
#' \item{r2_hf_res}{expected r2 for each randomly subsetted dataset}
#' \item{summary_r2_hf}{r2 mean and sd for each number of subsetted loci}
#' \item{nobs}{number of observations}
#' \item{nloc}{number of markers}
#' 
#' 
#' @references
#' Slate, J., David, P., Dodds, K. G., Veenvliet, B. A., Glass, B. C., Broad, T. E., & McEwan, J. C. (2004). 
#' Understanding the relationship between the inbreeding coefficient 
#' and multilocus heterozygosity: theoretical expectations and empirical data. Heredity, 93(3), 255-265.
#' 
#' Szulkin, M., Bierne, N., & David, P. (2010). HETEROZYGOSITY-FITNESS CORRELATIONS: A TIME FOR REAPPRAISAL. 
#' Evolution, 64(5), 1202-1217.
#' 
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) 
#'        
#' @examples
#' data(mouse_msats)
#' genotypes <- convert_raw(mouse_msats, miss = NA)
#' (out <- r2_hf(genotypes, subsets = c(2,4,6,8,10,12), nboot = 1000, type = "msats", 
#'               parallel = FALSE))
#' plot(out)
#' @export
#'
#'

r2_hf <- function(genotypes, subsets = NULL, nboot = 100, type = c("msats", "snps"), 
                  parallel = FALSE, ncores = NULL) {
    
#     if (!(steps > 1) | (steps > ncol(genotypes))) {
#         stop("steps have to be at least two and smaller or equal than the number of markers used")
#     }
    gtypes <- as.matrix(genotypes)
    
    # check for subset sequence
    if ((sum(subsets < 2)) != 0) stop("You cannot subset less than 2 markers")
    if ((sum(subsets > ncol(genotypes)))!= 0) stop("The number of subsetted markers cannot exceed the overall number of markers")
    if (any(subsets%%1 != 0)) stop("All subsets have to be specified by integers")
    
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
            # according to the miller paper, negative g2s are set to r2 = 0.
            if (g2 < 0) return(r2 <- 0)
            var_sh <- var(sMLH(gtypes))
            r2 <- g2/var_sh
            r2
        }
    }
    
    if (type == "snps") {
        calc_r2 <- function(gtypes) {
            g2 <- g2_snps(gtypes)[["g2"]]
            # according to the miller paper, negative g2s are set to r2 = 0.
            if (g2 < 0) return(r2 <- 0)
            var_sh <- var(sMLH(gtypes))
            r2 <- g2/var_sh
        }
    }
    
    # calculate r2 for full data
    r2_hf_full <- calc_r2(gtypes)
    
    # check if nboot = 0
    if ((nboot <= 0) | (is.null(subsets))) {
        res <- list(call = match.call(),
                    r2_hf_full = r2_hf_full,
                    r2_hf_res = NA,
                    summary_r2_hf = NA,
                    nobs = nrow(genotypes), 
                    nloc = ncol(genotypes))
        
        class(res) <- "inbreed"
        return(res)
    }
    
    # sorting
    subsets <- sort(subsets)
    
    # initialise
    all_r2 <- matrix(nrow = nboot, ncol = length(subsets))
    
    calc_r2_sub <- function(gtypes, i) {
        loci <- sample((1:ncol(gtypes)), i)
        out <- calc_r2(gtypes[, loci])
        out
    }
    
    # counter
    step_num <- 1
    
    if (parallel == FALSE) {
        
        for (i in subsets) {
            
            cat("\n", "Iterating subset number ", step_num, " from ", length(subsets), sep = "")
            if (step_num == length(subsets)) {
                cat("\n", "Last subset!", sep = "")
            }
            all_r2[, step_num] <- replicate(nboot, calc_r2_sub(gtypes, i))
            step_num <- step_num + 1
        }
        
    } else if (parallel == TRUE) {
        
        # define with counter for parallelized bootstraps
        calc_r2_sub_parallel <- function(boot_num, gtypes, i) {
            loci <- sample((1:ncol(gtypes)), i)
            out <- calc_r2(gtypes[, loci])
            out
        }
        
        for (i in subsets) {
            
            cat("\n", "Iterating subset number ", step_num, " from ", length(subsets), sep = "")
            if (step_num == length(subsets)) {
                cat("\n", "Last subset!", sep = "")
            }
            
            if (is.null(ncores)) {
                ncores <- parallel::detectCores()-1
                warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
            }
            
            cl <- parallel::makeCluster(ncores)
            all_r2[, step_num] <- parallel::parSapply(cl, 1:nboot, calc_r2_sub_parallel, gtypes, i)
            parallel::stopCluster(cl)

            step_num <- step_num + 1
        }
    }
    
    # expected r2 per subset
    r2_hf_res <- data.frame(r2 = c(all_r2), nloc = factor(rep(subsets, each = nboot)))
    
    # mean and sd per number of loci
    summary_r2_hf <- as.data.frame(as.list(aggregate(r2 ~ nloc, data = r2_hf_res, 
                                FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                                    sd = sd(x, na.rm = TRUE)))))
    names(summary_r2_hf) <- c("nloc", "Mean", "SD")
    
    res <- list(call = match.call(),
                r2_hf_full = r2_hf_full,
                r2_hf_res = r2_hf_res,
                summary_r2_hf = summary_r2_hf,
                nobs = nrow(genotypes), 
                nloc = ncol(genotypes))
    
    class(res) <- "inbreed"
    return(res)
    
}