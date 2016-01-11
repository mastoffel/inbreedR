#' Expected r2 between standardized multilocus heterozygosity (h) and inbreeding level (f)
#' 
#' 
#'
#' @param genotypes \code{data.frame} with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote), 1 (heterozygote) and \code{NA} (missing)
#' @param type specifies g2 formula to take. Type "snps" for large datasets and "msats" for smaller datasets.
#' @param nboot number of bootstraps over individuals to estimate a confidence interval
#'        around r2(h, f)
#' @param parallel Default is FALSE. If TRUE, bootstrapping and permutation tests are parallelized 
#' @param ncores Specify number of cores to use for parallelization. By default,
#'        all available cores but one are used.
#' @param CI confidence interval (default to 0.95)
#' @param subsets deprecated. a vector specifying the sizes of marker-subsets to draw. For a subset of 20 markers, subsets = c(2, 5, 10, 15, 20) could
#'        be a reasonable choice. The minimum subset size is 2 and the maximum is the number of markers in the data. 
#' @param nboot_loci deprecated. number of re-draws per subset of loci.
#'        
#' 
#' @return 
#' \item{call}{function call.}
#' \item{r2_hf_full}{expected r2 between inbreeding and sMLH for the full dataset}
#' \item{r2_hf_boot}{expected r2 values from bootstrapping over individuals}
#' \item{CI_boot}{confidence interval around the expected r2}
#' \item{r2_hf_res}{expected r2 for each randomly subsetted dataset}
#' \item{summary_r2_hf_res}{r2 mean and sd for each number of subsetted loci}
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
#' genotypes <- convert_raw(mouse_msats)
#' (out <- r2_hf(genotypes, nboot = 100, type = "msats", parallel = FALSE, nboot_loci = 100))
#' plot(out)
#' @export
#'
#'

r2_hf <- function(genotypes, type = c("msats", "snps"), nboot = NULL, 
                  parallel = FALSE, ncores = NULL, CI = 0.95, subsets = NULL, nboot_loci = 100) {
    
    
    # deprecate subsetting
    if (!is.null(subsets)) {
        warning("argument subsets is deprecated and will be deleted in the next package version. When you use it anyway, bear in mind
                that the variance of the estimates is biased due to subsetting from a finite number of markers.", 
                call. = FALSE)
    }
        
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
    
    # define g2 function
    if (type == "msats") {
        g2_fun <- g2_microsats
    } else if (type == "snps") {
        g2_fun <- g2_snps
    }
    
    # define calculation of expected r2
    calc_r2 <- function(gtypes) {
            g2 <- g2_fun(gtypes)[["g2"]]
            # according to the miller paper, negative g2s are set to r2 = 0.
            if (g2 < 0) return(r2 <- 0)
            var_sh <- stats::var(sMLH(gtypes))
            r2 <- g2/var_sh
            r2
    }
    
    # calculate r2 for full data
    r2_hf_full <- calc_r2(gtypes)
    
    # initialise r2_hf_boot
    r2_hf_boot <- NA
    CI_boot <- NA
    # bootstrapping over r2?
    if (!is.null(nboot)) {
        if (nboot < 2) stop("specify at least 2 bootstraps with nboot to 
                            estimate a confidence interval")
        # initialise
        r2_hf_boot <- matrix(nrow = nboot)
        
        # bootstrap function
        calc_r2_hf_boot <- function(boot, gtypes) {
            inds <- sample(1:nrow(gtypes), replace = TRUE)
            out <- calc_r2(gtypes[inds, ])
            # notifications
            if (boot %% 20 == 0) cat("\n", boot, "bootstraps over individuals done")
            if (boot == nboot) cat("\n", "### bootstrapping over individuals finished! ###")
            # result
            return(out)
        }
        
        # parallel bootstrapping ?
        if (parallel == FALSE) r2_hf_boot <- sapply(1:nboot, calc_r2_hf_boot, gtypes)
        
        if (parallel == TRUE) {
            if (is.null(ncores)) {
                ncores <- parallel::detectCores()-1
                warning("No core number specified: detectCores() is used to detect the number 
                        of \n cores on the local machine")
            }
            # run cluster
            cl <- parallel::makeCluster(ncores)
            parallel::clusterExport(cl, c("g2_microsats", "g2_snps", "sMLH"), envir = .GlobalEnv)
            r2_hf_boot <- parallel::parSapply(cl, 1:nboot, calc_r2_hf_boot, gtypes)
            parallel::stopCluster(cl)
        }
        r2_hf_boot <- c(r2_hf_full, r2_hf_boot)
        CI_boot <- stats::quantile(r2_hf_boot, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
    }
    
    # subsetting loci ?
    if ((nboot_loci <= 0) | (is.null(subsets))) {
        res <- list(call = match.call(),
                    r2_hf_full = r2_hf_full,
                    r2_hf_boot = r2_hf_boot,
                    CI_boot = CI_boot,
                    r2_hf_res = NA,
                    summary_r2_hf_res = NA,
                    nobs = nrow(genotypes), 
                    nloc = ncol(genotypes))
        
        class(res) <- "inbreed"
        return(res)
    }
    
    # sorting
    subsets <- sort(subsets)
    # initialise
    all_r2 <- matrix(nrow = nboot_loci, ncol = length(subsets))

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
            if (step_num == length(subsets)) cat("\n", "Last subset!", sep = "")
            
            all_r2[, step_num] <- replicate(nboot_loci, calc_r2_sub(gtypes, i))
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
            if (step_num == length(subsets)) cat("\n", "Last subset!", sep = "")
            
            if (is.null(ncores)) {
                ncores <- parallel::detectCores()-1
                warning("No core number specified: detectCores() is used to detect the number 
                        of \n cores on the local machine")
            }
            
            cl <- parallel::makeCluster(ncores)
            parallel::clusterExport(cl, c("g2_microsats", "g2_snps", "sMLH"), envir = .GlobalEnv)
            all_r2[, step_num] <- parallel::parSapply(cl, 1:nboot_loci, calc_r2_sub_parallel, gtypes, i)
            parallel::stopCluster(cl)

            step_num <- step_num + 1
        }
    }
    
    # expected r2 per subset
    r2_hf_res <- data.frame(r2 = c(all_r2), nloc = factor(rep(subsets, each = nboot_loci)))
    
    # mean and sd per number of loci
    summary_r2_hf_res <- as.data.frame(as.list(stats::aggregate(r2 ~ nloc, data = r2_hf_res, 
                                FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                                    sd = stats::sd(x, na.rm = TRUE),
                                                    CI_boot = stats::quantile(
                                                        x, c((1-CI)/2,1-(1-CI)/2), 
                                                        na.rm=TRUE)))))
    names(summary_r2_hf_res) <- c("nloc", "Mean", "SD", "CI_lower", "CI_upper")
    
    res <- list(call = match.call(),
                r2_hf_full = r2_hf_full,
                r2_hf_boot = r2_hf_boot,
                CI_boot = CI_boot,
                r2_hf_res = r2_hf_res,
                summary_r2_hf_res = summary_r2_hf_res,
                nobs = nrow(genotypes), 
                nloc = ncol(genotypes))
    
    class(res) <- "inbreed"
    return(res)
    
}