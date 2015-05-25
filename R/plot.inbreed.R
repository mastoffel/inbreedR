#' Plot an inbreed object
#' Plots the distribution of g2 estimates from bootstrapping.
#'
#' @param x An inbreed object.
#' @param main Plot title
#' @param xlab X-axis label
#' @param ylab y-axis label
#' @param \dots Additional arguments to the hist() function.
#'
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @seealso \link{g2_snps}, \link{g2_microsats}
#'
#'
#' @export
#'

plot.inbreed <- function(x, main = "g2 bootstrap values", xlab = "g2", ylab = "counts", ...) {
    
    if (!is.null(x$g2)) {
        boot_hist <- function(g2, g2_boot, CI.l, CI.u, ...) {
                # y position of confidence band
                v.pos <- max((hist(g2_boot, plot = FALSE, ...))$counts)
                # plot
                hist(g2_boot, ylim = c(0, v.pos*1.5), main = main, xlab = xlab, ylab = ylab, ...)
                lines(x = c(g2, g2), y = c(0, v.pos * 1.15), lwd = 2.5, col = "grey", lty = 5)
                arrows(CI.l, v.pos*1.15, CI.u, v.pos*1.15,
                       length=0.3, angle=90, code=3, lwd = 2.5, col = "black")
                points(g2, v.pos*1.15, cex = 1.2, pch = 19, col = "red")
                legend("topleft", pch = 19, cex = 1, bty = "n", col = c("red"),
                       c("g2 with CI"), box.lty = 0)
        }

        boot_hist(g2 = x$g2, g2_boot = x$g2_boot, CI.l = unname(x$CI_boot[1]), CI.u = unname(x$CI_boot[2]))
    }
    

      if (!is.null(x$sMLH_perm_mean)){
          perm_est <- mean(x$sMLH_perm_mean)
          perm_est_var <- var(x$sMLH_perm_var)
          range <- seq(min(x$sMLH_perm_mean), max(x$sMLH_perm_mean), length = 1000)
 
          hist(x$emp_sMLH, freq = FALSE)
          mean_het <- mean(x$emp_sMLH)
          sd_het <- sd(x$emp_sMLH)
          curve(dnorm(x, mean=mean_het, sd=sd_het), 
                col="darkblue", lwd=2, add=TRUE, yaxt="n")
      }
    
    g <- x$emp_sMLH
    h <- hist(g, breaks=14, density=100, col="lightgray", xlab="Accuracy", main="Overall") 
    xfit<-seq(min(g),max(g),length=10000) 
    yfit<-dnorm(xfit,mean=mean(x$sMLH_perm_mean),sd=mean(x$sMLH_perm_var)) 
    yfit <- yfit*diff(h$mids[1:2])*length(g) 
    lines(xfit, yfit, col="black", lwd=2)
    
}
