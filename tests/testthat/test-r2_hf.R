library(inbreedR)
context("r2_het_inbreeding")

# loading microsat data
data(mouse_msats)
msats <- convert_raw(mouse_msats, miss = NA)

# loading snp data
data(mouse_snps)
snps <- mouse_snps[1:500]

test_that("nboot = 0 works", {
    expect_equal(r2_hf(msats, subsets = c(2,4,10), nboot = 0, type = "msats")$r2_hf_full, 0.2797, tolerance = 0.001)
    expect_equal(r2_hf(snps, subsets = c(2,4,10), nboot = 0, type = "snps")$r2_hf_full, 0.6698, tolerance = 0.001)
})

test_that("subsets = NULL works", {
    expect_equal(r2_hf(msats, nboot = 0, type = "msats")$r2_hf_full, 0.2797, tolerance = 0.001)
    expect_equal(r2_hf(snps, nboot = 0, type = "snps")$r2_hf_full, 0.6698, tolerance = 0.001)
    expect_equal(r2_hf(msats, nboot = 5, type = "msats")$r2_hf_full, 0.2797, tolerance = 0.001)
    expect_equal(r2_hf(snps, nboot = 5, type = "snps")$r2_hf_full, 0.6698, tolerance = 0.001)
})

test_that("Matrix input works", {
    expect_equal(r2_hf(as.matrix(msats), subsets = c(2,4,10), nboot = 0, type = "msats")$r2_hf_full, 0.2797, tolerance = 0.001)
    expect_equal(r2_hf(as.matrix(snps), subsets = c(2,4,10), nboot = 0, type = "snps")$r2_hf_full, 0.6698, tolerance = 0.001)
})

