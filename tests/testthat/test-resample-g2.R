library(inbreedR)
context("resample_g2")

# loading microsat data
data(mouse_msats)
msats <- convert_raw(mouse_msats, miss = NA)

# loading snp data
data(mouse_snps)
snps <- mouse_snps[1:500]

test_that("g2 point estimate is correct", {
    expect_equal(round(resample_g2(msats, nboot = 0, subsets = NULL, type = "msats")$g2, 8), 0.02179805)
    expect_equal(round(resample_g2(snps, nboot = 0, subsets= NULL, type = "snps")$g2, 8), 0.02270884)
})

test_that("No bootstraps lead to NA´s in output", {
    expect_equal(is.na(resample_g2(msats, nboot = 0, subsets = NULL, type = "msats")$all_g2_res), TRUE)
    expect_equal(is.na(resample_g2(snps, nboot = 0, subsets = NULL, type = "snps")$all_g2_res), TRUE)
})

test_that("subset specifications work", {
    expect_equal(sum(is.na(resample_g2(msats, nboot = 5, 
                       subsets = c(2, 4, 6), type = "msats")$summary_all_g2$SD)), 0)
    expect_equal(sum(is.na(resample_g2(msats, nboot = 5, 
                                      subsets = c(2, 4, 6), type = "msats")$summary_all_g2$Mean)), 0)
})

test_that("subset exceeds marker number throws error", {

    expect_error(resample_g2(msats, nboot = 2, 
                            subsets = c(2, 4, 6, 1000), type = "msats"))
})

