#' Expected r2 between inbreeding level (f) and fitness (W)
#'
#' @param genotypes data.frame or matrix with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote) and 1 (heterozygote)
#' @param trait vector of any type which can be specified in R's glm() function
#' @param family distribution of the trait. Default is gaussian. For other distributions, just naming the distribution
#'        (e.g. binomial) will use the default link function (see ?family). Specifying another
#'        link function can be done in the same way as in the glm() function. A binomial distribution with 
#'        probit instead of logit link would be specified with family = binomial(link = "probit") 
#' @param subsets a vector specifying the sizes of subsets to draw. For a subset of 20 markers, subsets = c(2, 5, 10, 15, 20) could
#'        be a reasonable choice. The minimum subset size is 2 and the maximum is the number of markers in the data.
#' @param nboot number re-draws per subset size.
#' @param type specifies g2 formula to take. Type "snps" for large datasets and "msats" for smaller datasets.
#' 
#' @return 
#' \item{call}{function call.}
#' \item{exp_r2_full}{expected r2 between inbreeding and sMLH for the full dataset}
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
#' \dontrun{data(mouse_msats)
#' genotypes <- convert_raw(mouse_msats, miss = NA)
#' (out <- r2_Wf(genotypes, trait, family = gaussian, subsets = c(2,4,6,8,10,12), nboot = 100, type = "msats"))
#' plot(out)
#' }
#' @export
#'
#'

r2_Wf <- function(genotypes, trait, family = gaussian, subsets = NULL, nboot = 100, type = c("msats", "snps")) {
    
    # check for subset sequence
    if (!is.null(subsets)) {
        if ((sum(subsets < 2)) != 0) stop("You cannot subset less than 2 markers")
        if ((sum(subsets > ncol(genotypes)))!= 0) stop("The number of subsetted markers cannot exceed the overall number of markers")
        if (any(subsets%%1 != 0)) stop("All subsets have to be specified by integers")
    }
    
    # check g2 function argument
    if (length(type) == 2){
        type <- "msats"
    } else if (!((type == "msats")|(type == "snps"))){
        stop("type argument needs to be msats or snps")
    } 
    # check if trait is a vector
    if (!(is.atomic(trait) || is.list(trait))) {
        stop("trait has to be a vector")
    }
    
    # check for same number of individuals
    if (!(length(trait) == nrow(genotypes))) {
        stop("trait and genotypes have to contain the same number of individuals")
    }
    
    # check for data type of trait
    
    
    # genotypes matrix
    genotypes <- as.matrix(genotypes)
    
    # g2 function
    if (type == "msats") {
        g2_fun <- g2_microsats
    }
    
    if (type == "snps") {
        g2_fun <- g2_snps
    }
    
    # r2_Wf function
    calc_r2_Wf <- function(genotypes, trait) {
        # calculate sMLH
        het <- sMLH(genotypes)
        # Regression of trait on heterozygosity
        mod <- glm(trait ~ het, family = family)
        # beta coefficient
        beta_Wf <- coef(mod)[2]
        # R2 Wf
        R2 <- cor(trait,predict(mod))^2
        # g2
        g2 <- g2_fun(genotypes)[["g2"]]
        # according to the miller paper, negative g2´s are set to r2 = 0.
        if (g2 < 0) return( r2_Wf_res <- 0)
        # squared correlation between inbreeding and the fitness trait
        r2_Wf_res <- (R2 * var(het)) / (g2 * mean(het)^2)
        # r2_Hf_res <- R2 / r2_Wf_res
    }
    
    # r2_Wf for the full dataset 
    r2_Wf_full <- calc_r2_Wf(genotypes, trait)
    
    # check if nboot = 0
    if ((nboot <= 0) | (is.null(subsets))) {
        res <- list(call = match.call(),
                    r2_Wf_full = r2_Wf_full,
                    r2_Wf_res= NA,
                    summary_r2_Wf = NA,
                    nobs = nrow(genotypes), 
                    nloc = ncol(genotypes))
        
        class(res) <- "inbreed"
        return(res)
    }
    
    # sorting
    subsets <- sort(subsets)

    # initialise
    all_r2 <- matrix(nrow = nboot, ncol = length(subsets))
    
    calc_r2_Wf_sub <- function(genotypes, i) {
        loci <- sample((1:ncol(genotypes)), i)
        out <- calc_r2_Wf(genotypes[, loci], trait)
        out
    }
    
    # counter
    step_num <- 1
    for (i in subsets) {
        
        cat("\n", "Iterating subset number ", step_num, " from ", length(subsets), sep = "")
        if (step_num == length(subsets)) {
            cat("\n", "Last subset!", sep = "")
        }
        
        all_r2[, step_num] <- replicate(nboot, calc_r2_Wf_sub(genotypes, i))
        
        step_num <- step_num + 1
    }
    
    # expected r2 per subset
    r2_Wf_res <- data.frame(r2 = c(all_r2), nloc = factor(rep(subsets, each = nboot)))
    
    # mean and sd per number of loci
    summary_r2_Wf <- as.data.frame(as.list(aggregate(r2 ~ nloc, data = r2_Wf_res, 
                                                      FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                                                          sd = sd(x, na.rm = TRUE)))))
    names(summary_r2_Wf) <- c("nloc", "Mean", "SD")
    
    res <- list(call = match.call(),
                r2_Wf_full  = r2_Wf_full,
                r2_Wf_res = r2_Wf_res,
                summary_r2_Wf = summary_r2_Wf,
                nobs = nrow(genotypes), 
                nloc = ncol(genotypes))
    class(res) <- "inbreed"
    return(res)
    
}