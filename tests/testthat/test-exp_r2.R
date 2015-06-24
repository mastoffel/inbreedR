library(inbreedR)
context("expected r2")

# loading microsat data
data(mouse_msats)
msats <- convert_raw(mouse_msats, miss = NA)

# loading snp data
data(mouse_snps)
snps <- mouse_snps[1:500]

test_that("nboot = 0 works", {
    expect_equal(exp_r2(msats, subsets = c(2,4,10), nboot = 0, type = "msats")$exp_r2_full, 0.2797, tolerance = 0.001)
    expect_equal(exp_r2(snps, subsets = c(2,4,10), nboot = 0, type = "snps")$exp_r2_full, 0.6698, tolerance = 0.001)
})