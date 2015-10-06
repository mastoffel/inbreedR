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
    
    
    # check if its r2_hf
    if(!is.null(x$r2_hf_res)) {
        cat("\n\n", "Calculation of expected r2 between inbreeding level (f) and heterozygosity (sMLH)", "\n",
                    "---------------------------------------------------------------------------------", "\n\n", sep = "")
        cat("\n", "Data: ", x$nobs, " observations at ", x$nloc, " markers", "\n",
            "Function call = ", deparse(x$call), "\n\n",
            sep = "")
        cat("Expected r2 based on all markers: ", x$r2_hf_full, "\n\n")
        cat("Confidence interval for r2 based on bootstrapping over individuals:", "\n")
        print(x$CI_boot)
        cat("\n")
        if(is.data.frame(x$summary_r2_hf_res)){
            cat("Average expected r2 of each marker subset: ", "\n\n")
            print(format(x$summary_r2_hf_res, digits = 2))
        }
    }
    
    # check if its r2_Wf
        if(!is.null(x$r2_Wf_full)) {
            cat("\n\n", "Calculation of expected r2 between inbreeding level (f) and fitness (W)", "\n",
                "---------------------------------------------------------------------------------", "\n\n", sep = "")
            cat("\n", "Data: ", x$nobs, " observations at ", x$nloc, " markers", "\n",
                "Function call = ", deparse(x$call), "\n\n",
                sep = "")
            cat("Expected r2 based on all markers: ", x$r2_Wf_full, "\n\n")
            cat("Confidence interval for r2 based on bootstrapping over individuals:", "\n")
            print(x$CI_boot)
            cat("\n")
        }
    
    
    
    # no subsetting yet
#     if(!is.null(x$r2_Wf_res)) {
#         cat("\n\n", "Calculation of expected r2 between inbreeding level (f) and fitness (W)", "\n",
#             "---------------------------------------------------------------------------------", "\n\n", sep = "")
#         cat("\n", "Data: ", x$nobs, " observations at ", x$nloc, " markers", "\n",
#             "Function call = ", deparse(x$call), "\n\n",
#             sep = "")
#         cat("Expected r2 based on all markers: ", x$r2_Wf_full, "\n\n")
#         if(is.data.frame(x$summary_r2_Wf)){
#             cat("Average expected r2 of each marker subset: ", "\n\n")
#             print(format(x$summary_r2_Wf, digits = 2))
#         }
#     }
    
    # check its resample_g2
    if(!is.null(x$all_g2_res)) {
        cat("\n\n", "Calculation of identity disequilibrium (g2) for an increasing number of markers", "\n",
            "---------------------------------------------------------------------------------", "\n\n", sep = "")
        cat("\n", "Data: ", x$nobs, " observations at ", x$nloc, " markers", "\n",
            "Function call = ", deparse(x$call), "\n\n",
            sep = "")
        cat("g2 estimate based on all markers: ", x$g2_full, "\n\n")
        if(is.data.frame(x$summary_all_g2)){
            cat("Average g2 of each marker subset: ", "\n\n")
            print(format(x$summary_all_g2, digits = 2))
        }
    }

    # check if HHC
    if(!is.null(x$HHC_vals)) {
        cat("\n\n", "heterozygosity-heterozygosity correlations", "\n",
                    "------------------------------------------", "\n", sep = "")
        cat("\n", "Data: ", x$nobs, " observations at ", x$nloc, " markers", "\n",
            "Function call = ", deparse(x$call), "\n",
            sep = "")
        cat("\n", "HHC Mean : " ,round(mean(x$HHC_vals, na.rm = TRUE), 3),
            "\n", "HHC SD: " , round(stats::sd(x$HHC_vals, na.rm = TRUE), 3), 
            "\n", "HHC CI: [", round(x$CI_HHC[1],3), ", ", round(x$CI_HHC[2],3), "]", sep = "")
    }
}
