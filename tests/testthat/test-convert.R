library(inbreedR)
context("convert_raw")

# load data with NA as missing
data(mouse_msats)
data <- mouse_msats

## constructing different data sets

# data with NA missing
data_NA <- data
# data with -1 missing
data[is.na(data)] <- -1
data_NUM <- data

# character data frame (for SNP data)
bp <- c("A","T","G","C")
m <- as.data.frame(matrix(data = bp, nrow = 100, ncol = 100))
m[] <- lapply(m, as.character)
m[] <- lapply(m, sample)
m[sample(1:100, 10), sample(1:100, 10)] <- NA
data_CHAR <- m

# convert raw
# different missing data
data_conv_NA <- convert_raw(data_NA, miss = NA)
data_conv_NUM <- convert_raw(data_NUM, miss = -1)
# Characters as genotypes
data_conv_CHAR <- convert_raw(data_CHAR, miss = NA)

test_that("Converted data contains half the number columns", {
    expect_equal(ncol(convert_raw(data_NA, miss = NA)), ncol(data_NA)/2)
    expect_equal(ncol(convert_raw(data_NUM, miss = -1)), ncol(data_NUM)/2)
    expect_equal(ncol(convert_raw(data_CHAR, miss = NA)), ncol(data_CHAR)/2)
})

test_that("Converted data contains just 1, 0 and -1", {
    expect_equal(sum(!((data_conv_NA == 1) | (data_conv_NA == 0) | (is.na(data_conv_NA)))), 0)
    expect_equal(sum(!((data_conv_NUM == 1) | (data_conv_NUM == 0) | (is.na(data_conv_NUM)))), 0)
    expect_equal(sum(!((data_conv_CHAR == 1) | (data_conv_CHAR == 0) | (is.na(data_conv_CHAR)))), 0)
})






# # construct datasets
# library(dplyr)
# raw_data <- sample(c(1,0), 10000, replace = TRUE, prob = c(0.4,0.7))
# genotypes <- as.data.frame(matrix(raw_data, ncol = 1000))
# # genotypes_alleles <- as.data.frame(matrix(data = NA, nrow = nrow(genotypes), ncol = ncol(genotypes) * 2))
# 
# genotypes_alleles <- as.data.frame(matrix(data = NA, ncol = ncol(genotypes) * 2, nrow = nrow(genotypes)))
# second_col <- vector(1:ncol(genotypes_alleles), by = 2)
# 
# 
# make_alleles <- function(locus) {
#     locus <- as.matrix(locus)
#     alleles <- matrix(rep(locus, 2), ncol = 2)
#     # create heterozygotes
#     alleles[alleles[,1] == 1, 1] <- 2
#     return(as.data.frame(alleles))
# }
# 
# out <- apply(genotypes, 2, make_alleles)
# genotypes_alleles <- do.call(cbind, t(apply(genotypes, 2, make_alleles)))
# 
# genotypes_alleles <- matrix(unlist(lapply(genotypes, make_alleles)), ncol = ncol(genotypes) * 2)
# 
# out <- convert_raw(genotypes_alleles , miss = NA)

