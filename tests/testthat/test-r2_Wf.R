context("r2_Wf")

# loading microsat data
data(mouse_msats)
msats <- convert_raw(mouse_msats)

# loading snp data
data(mouse_snps)
snps <- mouse_snps[1:500]

# loading bodyweight
data("bodyweight")

test_that("point estimate correct", {
    expect_equal(r2_Wf(msats, bodyweight, family = gaussian)$r2_Wf_full, 0.433, tolerance = 0.001)
    expect_equal(r2_Wf(snps, bodyweight, family = gaussian, type = "snps")$r2_Wf_full, 0.1287353, tolerance = 0.001)
})

test_that("matrix input works", {
    expect_equal(r2_Wf(as.matrix(msats), bodyweight, family = gaussian)$r2_Wf_full, 0.433, tolerance = 0.001)
    expect_equal(r2_Wf(as.matrix(snps), bodyweight, family = gaussian, type = "snps")$r2_Wf_full, 0.1287353, tolerance = 0.001)
})
