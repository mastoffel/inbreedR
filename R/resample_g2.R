#' Identity disequilibrium (g2) for different marker subsets
#' 
#' 
#' @param genotypes data.frame with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote) and 1 (heterozygote)
#' @param subsets a vector specifying the sizes of subsets to draw. For a subset of 20 markers, subsets = c(2, 5, 10, 15, 20) could
#'        be a reasonable choice. The minimum subset size is 2 and the maximum is the number of markers in the data.
#' @param nboot number re-draws per subset size.
#' @param type specifies g2 formula to take. Type "snps" for large datasets and "msats" for smaller datasets.
#' 
#' @return 
#' \item{call}{function call.}
#' \item{all_g2_res}{vector of g2 values for each randomly subsetted dataset}
#' \item{summary_exp_r2}{g2 mean and sd for each number of subsetted loci}
#' \item{nobs}{number of observations}
#' \item{nloc}{number of markers}
#' 
#' @references
#' Hoffman, J.I., Simpson, F., David, P., Rijks, J.M., Kuiken, T., Thorne, M.A.S., Lacey, R.C. & Dasmahapatra, K.K. (2014) High-throughput sequencing reveals inbreeding depression in a natural population.
#' Proceedings of the National Academy of Sciences of the United States of America, 111: 3775-3780. Doi: 10.1073/pnas.1318945111
#'
#' 
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) 
#'        
#' @examples
#' data(mouse_msats)
#' genotypes <- convert_raw(mouse_msats, miss = NA)
#' (out <- resample_g2(genotypes, subsets = c(2,4,6,8,10,12), nboot = 1000, type = "msats"))
#' plot(out)
#' @export
#'
#'

resample_g2 <- function(genotypes, subsets = NULL, nboot = 100, type = c("msats", "snps")) {
    
    genotypes <- as.matrix(genotypes)
    
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
    
    # assign g2 function
    if (type == "msats"){
        g2_fun <- g2_microsats
    } else {
        g2_fun <- g2_snps
    }
    
    # initialise
    nloc <- ncol(genotypes)
    all_g2 <- matrix(data = NA, nrow = nboot, ncol = length(subsets))
    
    # subsampling and calculating g2
    sample_genotypes <- function(genotypes, num_subsamp) {
        ind <- sample(1:ncol(genotypes), num_subsamp)
        g2 <- g2_fun(genotypes[, ind])[["g2"]]
        g2
    }
    
    step_num <- 1
    
    for (i in subsets) {
        
        cat("\n", "Iterating subset number ", step_num, " from ", length(subsets), sep = "")
        if (step_num == length(subsets)) {
            cat("\n", "Last subset!", sep = "")
        }
        
        all_g2[, step_num] <- replicate(nboot, sample_genotypes(genotypes, i))
      
        step_num <- step_num + 1
    }
    
    # variable names are number of markers used
    #all_g2 <- as.data.frame(all_g2)
    
    # expected r2 per subset
    all_g2_res <- data.frame(g2 = c(all_g2), nloc = factor(rep(subsets, each = nboot)))
    
    # mean and sd per number of loci
    summary_all_g2 <- as.data.frame(as.list(aggregate(g2 ~ nloc, data = all_g2_res, 
                                FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                                    sd = sd(x, na.rm = TRUE)))))
    
    res <- list(call = match.call(),
                all_g2_res = all_g2_res,
                summary_all_g2 = summary_all_g2,
                nobs = nrow(genotypes), 
                nloc = ncol(genotypes))
    
    class(res) <- "inbreed"
    return(res)
    
}


