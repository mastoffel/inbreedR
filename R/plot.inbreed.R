#' Plot an inbreed object
#' 
#'
#' @param x An inbreed object.
#' @param plottype "boxplot" or "histogram" to plot the output of r2_hf() and to show either
#'        the boxplots through resampling of loci or the histogram from the bootstrapping of 
#'        r2 over individuals.
#' @param \dots Additional arguments to the hist() function for the g2 and the HHC functions. 
#'              Additional arguments to the boxplot() function for plotting the result of the r2_hf() function.
#'
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @seealso \link{g2_snps}, \link{g2_microsats}
#'
#' @export
#'

plot.inbreed <- function(x, plottype = c("boxplot", "histogram"), ...) {
    # check if its a g2 calculater
    if (!is.null(x[["g2"]])) {
        if (is.na(x$g2_se)) stop("Number of bootstraps specified in g2 function was 0, so there is nothing to plot")
        # save ellipsis args
        dots <- list(...)
        # make empty arguments list
        args1 <- list()
        if (!("main" %in% dots)) args1$main <- "g2 bootstrap distribution"
        if (!("xlab" %in% dots)) args1$xlab <- "g2"
        if (!("ylab" %in% dots)) args1$ylab <- "counts"
        # add names (will be argument names) to args1 values
        args1[names(dots)] <- dots
        
        boot_hist <- function(g2, g2_boot, CI.l, CI.u, args1) {
            
                # y position of confidence band
                v.pos <- max(do.call(graphics::hist, (c(list(x = g2_boot, plot = FALSE, warn.unused = FALSE), args1)))$counts)
                # plot
                do.call(graphics::hist, (c(list(x = g2_boot, ylim = c(0, v.pos*1.5)), args1))) 
                graphics::lines(x = c(g2, g2), y = c(0, v.pos * 1.15), lwd = 2.5, col = "black", lty = 5)
                graphics::arrows(CI.l, v.pos*1.15, CI.u, v.pos*1.15,
                       length=0.1, angle=90, code=3, lwd = 2.5, col = "black")
                graphics::points(g2, v.pos*1.15, cex = 1.2, pch = 19, col = "black")
                graphics::legend(x = "topleft", inset = 0.01, pch = 19, cex = 1, bty = "n", col = c("black"),
                       c("g2 with CI"), box.lty = 0)
        }
        boot_hist(g2 = x$g2, g2_boot = x$g2_boot, CI.l = unname(x$CI_boot[1]), 
                  CI.u = unname(x$CI_boot[2]), args1 = args1)
    }
    # check if its r2_hf
    if (as.character(x$call[[1]]) == "r2_hf") {
        # save ellipsis args
        dots <- list(...)
        # if plottype argument not specified, assign boxplot
        if (length(plottype) == 2) plottype <- plottype[1]

        # plotting resampling boxplots
        if (plottype == "boxplot") {
            
        if(!is.null(x$r2_hf_res)) {
        # check whether values available
        if(!is.data.frame(x$summary_r2_hf_res)) stop("No resampling done, so nothing to plot")
        
        # make empty arguments list
        args1 <- list()
        
        if (!("main" %in% dots)) args1$main <- "Expected r2 between f and sMLH"
        if (!("xlab" %in% dots)) args1$xlab <- "number of loci"
        if (!("ylab" %in% dots)) args1$ylab <- "r2"
        
        # add names (will be argument names) to args1 values
        args1[names(dots)] <- dots
        
        do.call(graphics::boxplot, c(list(formula = r2 ~ nloc, data = x$r2_hf_res,
                pch = 16, outcol ="black"), args1))
        }
        } else if (plottype == "histogram") {
            if (length(x$r2_hf_boot) < 2) stop("No bootstrapping done, so nothing to plot")
            
            # make empty arguments list
            args1 <- list()
            if (!("main" %in% dots)) args1$main <- "r2 bootstrap distribution"
            if (!("xlab" %in% dots)) args1$xlab <- "r2 (heterozygosity, inbreeding)"
            if (!("ylab" %in% dots)) args1$ylab <- "counts"
            # add names (will be argument names) to args1 values
            args1[names(dots)] <- dots
            
            boot_hist_r2_hf <- function(r2_hf_full, r2_hf_boot, CI.l, CI.u, args1) {
                
                # y position of confidence band
                v.pos <- max(do.call(graphics::hist, (c(list(x = r2_hf_boot, plot = FALSE, 
                                                             warn.unused = FALSE), args1)))$counts)
                # plot
                do.call(graphics::hist, (c(list(x = r2_hf_boot, ylim = c(0, v.pos*1.5)), args1))) 
                graphics::lines(x = c(r2_hf_full, r2_hf_full), y = c(0, v.pos * 1.15), lwd = 2.5, 
                                col = "black", lty = 5)
                graphics::arrows(CI.l, v.pos*1.15, CI.u, v.pos*1.15,
                                 length=0.1, angle=90, code=3, lwd = 2.5, col = "black")
                graphics::points(r2_hf_full, v.pos*1.15, cex = 1.2, pch = 19, col = "black")
                graphics::legend(x = "topleft", inset = 0.01, pch = 19, cex = 1, bty = "n", 
                                 col = c("black"), c("r2 with CI"), box.lty = 0)
            }
            boot_hist_r2_hf(r2_hf_full = x$r2_hf_full, r2_hf_boot = x$r2_hf_boot, 
                      CI.l = unname(x$CI_boot[1]), 
                      CI.u = unname(x$CI_boot[2]), args1 = args1)
        }
    }
    
    
    # check if its resample_g2
    if(!is.null(x$all_g2_res)) {
        # save ellipsis args
        dots <- list(...)
        # make empty arguments list
        args1 <- list()
        
        if (!("main" %in% dots)) args1$main <- "g2 for different marker subsets"
        if (!("xlab" %in% dots)) args1$xlab <- "number of loci"
        if (!("ylab" %in% dots)) args1$ylab <- "g2"
        
        # add names (will be argument names) to args1 values
        args1[names(dots)] <- dots
        
        do.call(graphics::boxplot, c(list(formula = g2 ~ nloc, data = x$all_g2_res,
                                pch = 16, outcol ="black"), args1))
    }
    
    # check if HHC
    if(!is.null(x$HHC_vals)) {
        # save ellipsis args
        dots <- list(...)
        # make empty arguments list
        args1 <- list()
        
        if (!("main" %in% dots)) args1$main <-  "heterozygosity-heterozygosity correlation distribution"
        if (!("xlab" %in% dots)) args1$xlab <- "correlation coefficient r"
        if (!("ylab" %in% dots)) args1$ylab <- "counts"
        
        # add names (will be argument names) to args1 values
        args1[names(dots)] <- dots
        
        boot_hist_HHC <- function(g2, g2_boot, CI.l, CI.u, args1) {
            # y position of confidence band
            v.pos <- max(do.call(graphics::hist, (c(list(x = g2_boot, plot = FALSE, warn.unused = FALSE), args1)))$counts)
            # plot
            do.call(graphics::hist, (c(list(x = g2_boot, ylim = c(0, v.pos*1.5)), args1))) 
            graphics::lines(x = c(g2, g2), y = c(0, v.pos * 1.15), lwd = 2.5, col = "black", lty = 5)
            graphics::arrows(CI.l, v.pos*1.15, CI.u, v.pos*1.15,
                   length=0.1, angle=90, code=3, lwd = 2.5, col = "black")
            graphics::points(g2, v.pos*1.15, cex = 1.2, pch = 19, col = "black")
            graphics::legend(x = "topleft", inset = 0.01, pch = 19, cex = 1, bty = "n", col = c("black"),
                   c("mean HHC with CI"), box.lty = 0)
        }
        boot_hist_HHC(g2 = mean(x$HHC_vals, na.rm = TRUE), g2_boot = x$HHC_vals, 
                  CI.l = x$CI_HHC[1], CI.u = x$CI_HHC[2], args1 = args1)
    }
        
}
