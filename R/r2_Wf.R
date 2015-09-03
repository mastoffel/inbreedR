#' Expected r2 between inbreeding level (f) and fitness (W)
#'
#' @param genotypes data.frame with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote), 1 (heterozygote) and NA (missing)
#' @param trait vector of any type which can be specified in R's glm() function
#' @param family distribution of the trait. Default is gaussian. For other distributions, just naming the distribution
#'        (e.g. binomial) will use the default link function (see ?family). Specifying another
#'        link function can be done in the same way as in the glm() function. A binomial distribution with 
#'        probit instead of logit link would be specified with family = binomial(link = "probit") 
#' @param type specifies g2 formula to take. Type "snps" for large datasets and "msats" for smaller datasets.
#' 
#' @return 
#' \item{call}{function call.}
#' \item{exp_r2_full}{expected r2 between inbreeding and sMLH for the full dataset}
#' \item{nobs}{number of observations}
#' \item{nloc}{number of markers}
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
#' data(bodyweight)
#' genotypes <- convert_raw(mouse_msats)
#' 
#' (out <- r2_Wf(genotypes, bodyweight, family = "gaussian", type = "msats"))
#' 
#' 
#' @export
#'
#'

r2_Wf <- function(genotypes, trait, family = "gaussian", type = c("msats", "snps")) {
    
    # genotypes matrix
    genotypes <- as.matrix(genotypes)
    
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
        mod <- stats::glm(trait ~ het, family = family)
        # beta coefficient
        beta_Wh <- stats::coef(mod)[2]
        # R2 Wh
        R2 <- stats::cor(trait,stats::predict(mod))^2
        # g2
        g2 <- g2_fun(genotypes)[["g2"]]
        # according to the miller paper, negative g2s are set to r2 = 0.
        # if (g2 < 0) return( r2_Wf_res <- 0)
        # squared correlation between inbreeding and the fitness trait
        r2_Wf_res <- (R2 * stats::var(het, na.rm = TRUE)) / (g2 * mean(het, na.rm = TRUE)^2)
        # r2_Hf_res <- R2 / r2_Wf_res
    }
    
    # r2_Wf for the full dataset 
    r2_Wf_full <- calc_r2_Wf(genotypes, trait)
    
    res <- list(call = match.call(),
                r2_Wf_full  = r2_Wf_full,
                nobs = nrow(genotypes), 
                nloc = ncol(genotypes))
    class(res) <- "inbreed"
    return(res)
    
}