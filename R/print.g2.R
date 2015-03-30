#' Print a rpt object
#'
#' Displays the results a rpt object (i.e. the result of a rpt function call) an a nice form.
#'
#' @param x A g2 object return from one of the g2 functions.
#' @param \dots Additional arguments; none are used in this method.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @seealso \link{g2_snps}, \link{g2_microsats}, \link{plot}
#'
#' @export



print.g2 <- function(x, ...) {

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
