library(inbreedR)
context("g2 functions")

# loading microsat data
data(seal_microsats)
msats <- convert_raw(seal_microsats, miss = NA)

# loading snp data
data(mice_snp_genotypes)
snps <- mice_snp_genotypes[1:500]

test_that("g2 point estimates are correct", {
    expect_equal(round(g2_microsats(msats, nperm = 0, nboot = 0, CI = 0.95)$g2, 8), 0.00241214)
    expect_equal(round(g2_snps(snps, nperm = 0, nboot = 0, CI = 0.95)$g2, 8), 0.02270884)
    })

test_that("bootstrapping worked", {
    expect_equal(sum(is.na(g2_microsats(msats, nperm = 0, nboot = 20, CI = 0.95)$CI_boot)), 0) # CI calculated
    expect_equal(sum(is.na(g2_snps(snps, nperm = 0, nboot = 20, CI = 0.95)$CI_boot)), 0)
})

test_that("permutation worked", {
    expect_equal(sum(is.na(g2_microsats(msats, nperm = 20, nboot = 0, CI = 0.95)$CI_boot)), 0) # CI calculated
    expect_equal(sum(is.na(g2_microsats(msats, nperm = 20, nboot = 20, CI = 0.95)$CI_boot)), 0)
})