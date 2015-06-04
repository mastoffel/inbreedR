#' inbreedR: Workflows for analysing variance in inbreeding based on SNP or microsatellite markers.
#'
#' @description 
#' 
#' inbreedR contains the following functions:
#'
#' \link{g2_microsats}
#' \link{g2_snps}
#' \link{convert_raw}
#' \link{exp_r2}
#' \link{HHC}
#' \link{sMLH}
#' \link{plot.inbreed}
#'
#' @details 
#' 
#' It has three main goals:
#' 
#' \itemize{
#' \item Assessing identity disequilibria and the potential to detect heterozygosity-fitness correlations based on global effects
#' \item Providing insights on the sensitivity of these measures based on the number of molecular markers used
#' \item Implementing computationally efficient functions in a flexible environment for analysing 
#' inbreeding and HFC`s with both small and large datasets.
#' }
#' 
#' For a short introduction to inbreedR start with the vignette:
#' \code{browseVignettes(package = "inbreedR")}
#'
#' @author  Martin Stoffel (martin.adam.stoffel@@gmail.com), Mareike Esser (messer@@uni-bielefeld.de)
#'
#' @references
#' David, P., Pujol, B., Viard, F., Castella, V. and Goudet, J. (2007),
#' Reliable selfing rate estimates from imperfect population genetic data. Molecular Ecology,
#' 16: 2474
#'
#' Hoffman, J.I., Simpson, F., David, P., Rijks, J.M., Kuiken, T., Thorne, M.A.S., Lacey, R.C. & Dasmahapatra, K.K. (2014) High-throughput sequencing reveals inbreeding depression in a natural population.
#' Proceedings of the National Academy of Sciences of the United States of America, 111: 3775-3780. 
#'
#' @docType package
#' @name inbreedR
NULL
