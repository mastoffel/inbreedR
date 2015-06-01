#' Plot an inbreed object
#' 
#'
#' @param x An inbreed object.
#' @param \dots Additional arguments to the hist() function.
#'
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @seealso \link{g2_snps}, \link{g2_microsats}
#'
#' @import ggplot2
#' @export
#'

plot.inbreed <- function(x, ...) {
    # check if its a g2 calculater
    if (!is.null(x$g2)) {
        if (is.na(x$g2_se)) stop("Number of bootstraps specified in g2 function was 0, so there is nothing to plot")
        if (!(hasArg("main"))) main <- "g2 bootstrap values"
        if (!(hasArg("xlab")))   xlab <- "g2"
        if (!(hasArg("ylab"))) ylab <- "counts"
        
        boot_hist <- function(g2, g2_boot, CI.l, CI.u, ...) {
                # y position of confidence band
                v.pos <- max((hist(g2_boot, plot = FALSE, ...))$counts)
                # plot
                hist(g2_boot, ylim = c(0, v.pos*1.5), main = main, xlab = xlab, ylab = ylab, ...)
                lines(x = c(g2, g2), y = c(0, v.pos * 1.15), lwd = 2.5, col = "grey", lty = 5)
                arrows(CI.l, v.pos*1.15, CI.u, v.pos*1.15,
                       length=0.3, angle=90, code=3, lwd = 2.5, col = "black")
                points(g2, v.pos*1.15, cex = 1.2, pch = 19, col = "red")
                legend(x = "topleft", inset = 0.01, pch = 19, cex = 1, bty = "n", col = c("red"),
                       c("g2 with CI"), box.lty = 0)
        }
        boot_hist(g2 = x$g2, g2_boot = x$g2_boot, CI.l = unname(x$CI_boot[1]), CI.u = unname(x$CI_boot[2]))
    }
    # check if its exp_r2
    if(!is.null(x$exp_r2_res)) {
        ggplot(x$exp_r2_res, ggplot2::aes_string("nloc", "r2")) +
            geom_boxplot() +
            theme_bw(base_size = 16) +
            theme(plot.title=element_text(size=15, vjust=2)) +
            labs(title = "Expected r2 between f and sMLH")+
            xlab("number of loci") 
    }
    
    # check if HHC
    if(!is.null(x$HHC_vals)) {
        if (!(hasArg("main"))) main <- "heterozygosity-heterozygosity correlations"
        if (!(hasArg("xlab")))   xlab <- "correlation coefficient r"
        if (!(hasArg("ylab"))) ylab <- "counts"
        
        boot_hist <- function(g2, g2_boot, CI.l, CI.u, ...) {
            # y position of confidence band
            v.pos <- max((hist(g2_boot, plot = FALSE, ...))$counts)
            # plot
            hist(g2_boot, ylim = c(0, v.pos*1.5), main = main, xlab = xlab, ylab = ylab, ...)
            lines(x = c(g2, g2), y = c(0, v.pos * 1.15), lwd = 2.5, col = "grey", lty = 5)
            arrows(CI.l, v.pos*1.15, CI.u, v.pos*1.15,
                   length=0.3, angle=90, code=3, lwd = 2.5, col = "black")
            points(g2, v.pos*1.15, cex = 1.2, pch = 19, col = "red")
            legend(x = "topleft", inset = 0.01, pch = 19, cex = 1, bty = "n", col = c("red"),
                   c("mean HHC with CI"), box.lty = 0)
        }
        boot_hist(g2 = mean(x$HHC_vals, na.rm = TRUE), g2_boot = x$HHC_vals, CI.l = x$CI_HHC[1], CI.u = x$CI_HHC[2])
    }
        
}
