#' Print an inbreed object
#'
#' Displays the results a inbreed object.
#'
#' @param x An inbreed object from one of the inbreedR functions.
#' @param \dots Additional arguments; none are used in this method.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @seealso \link{g2_snps}, \link{g2_microsats}, \link{plot}
#'
#' @export



print.inbreed <- function(x, ...) {
    # check if its a g2 calculater
    if (!is.null(x$g2)) {
        if (as.character(x$call)[1] == "g2_microsats") {
                cat("\n\n", "Calculation of identity disequilibrium with g2 for microsatellite data", "\n",
                          "----------------------------------------------------------------------", "\n", sep = "")
        } else if (as.character(x$call)[1] == "g2_snps") {
                cat("\n\n", "Calculation of identity disequilibrium with g2 for SNP data", "\n",
                            "-----------------------------------------------------------", "\n", sep = "")
        }

        cat("\n", "Data: ", x$nobs, " observations at ", x$nloc, " markers", "\n",
            "Function call = ", deparse(x$call), "\n\n",
            sep = "")
        cat("g2 = ", x$g2, ", se = ", x$g2_se, "\n\n", sep  ="")
        cat("confidence interval", "\n")
        print(x$CI_boot)
        cat("\n")
        cat("p (g2 > 0) = ",  x$p_val, " (based on ", length(x$g2_permut), " permutations)",  sep = "")
    }
    
#     # check if its the variance comparer
#     if (!is.null(x$sMLH_perm_mean)){
#         df <- data.frame("mean" = c(mean(x$emp_sMLH), mean(x$sMLH_perm_mean), mean(x$sMLH_perm_var)),
#                "variance" = c(var(x$emp_sMLH), var(x$sMLH_perm_mean), var(x$sMLH_perm_var)))
#         row.names(df) <- c("sample", "perm_means", "perm_vars")
#         
#         cat("\n\n", "Comparison of empirical sMLH to a theoretical distribution in the absence of inbreeding:", sep = "")
#         cat("\n\n")
#         print(format(df, digits = 3, width = 6, scientific = FALSE))  
#     }
    
    
    # check if its exp_r2
    if(!is.null(x$exp_r2_res)) {
        cat("\n\n", "Calculation of expected r2 between f and sMLH for an increasing number of markers", "\n",
                    "---------------------------------------------------------------------------------", "\n\n", sep = "")
        cat("\n", "Data: ", x$nobs, " observations at ", x$nloc, " markers", "\n",
            "Function call = ", deparse(x$call), "\n\n",
            sep = "")
        print(format(x$summary_exp_r2, digits = 2))
    }

    
}
