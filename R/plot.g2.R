#' Plot a g2 object
#' Plots the distribution of g2 estimates from bootstrapping.
#'
#' @param x A g2 object.
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

plot.g2 <- function(x, main = "g2 bootstrap estimates", xlab = "g2", ylab = "counts", ...) {

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
