#' Checks the data for consistency with the inbreedR working format.
#' 
#' The inbreedR working format is an i * l genotype matrix, whereby each individual is a row 
#' and each column is a locus.
#' Heterozygosity at a given locus should be coded as 1, homozygosity as 0 and missing values
#' should be coded as NA.
#'
#' @param genotypes data.frame with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote), 1 (heterozygote) and NA (missing)
#' @param num_ind Number of individuals
#' @param num_loci Number of loci / markers
#'
#' @details Checks that (1) the genotype data just contains 3 elements, which is 0 for homozygote,
#'          1 for heterozygote and NA for missing data, (2) the number 
#'          of individuals corresponds to the number of rows and the number of loci corresponds to the
#'          number of columns, (3) the data type is numeric.
#'      .
#'          
#' @return
#' TRUE if the data format is correct, error message if any test failed
#'
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) 
#'   
#' @examples
#' data(mouse_msats)
#' # tranform raw genotypes into 0/1 format
#' genotypes <- convert_raw(mouse_msats)
#' # check data
#' check_data(genotypes, num_ind = 36, num_loci = 12)
#' 
#'
#' @export
#'


check_data <- function(genotypes, num_ind = NULL, num_loci = NULL) {

genotypes <- data.table::as.data.table(genotypes)
vals <- unique(unlist(lapply(genotypes, unique)))

if (length(unique(unlist(vals))) > 3) {
    stop("The data contains more than 3 elements (1, 0, missing)")
}

if (!(1 %in% (unique(unlist(vals))))) {
    stop("Heterozygosity at a locus is not coded as 1")
}

if (!(0 %in% (unique(unlist(vals))))) {
    stop("Homozygosity at a locus is not coded as 0")
}

if (!(NA %in% (unique(unlist(vals))))) {
    stop("Missing values have to be coded as NA")
}

if (!is.null(num_ind)) {
    if (num_ind != nrow(genotypes)) {
        stop("Number of rows is unequal to the number of individuals")
    }
}

if (!is.null(num_loci)) {
    if (num_loci != ncol(genotypes)) {
        stop("Number of columns is unequal to the number of loci")
    }
}

### check here for data types (todo) #### 
# classes <- apply(genotypes, 2, class)
# 
# if(any(classes != "numeric") & any(classes != "integer")) {
#     stop("At least some elements in the data are not numeric, 
#          check that you coded missing values with either NA or a number(except 0 and 1)")
# }

return(TRUE)

}
