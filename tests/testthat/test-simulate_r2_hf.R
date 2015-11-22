context("simulate_r2_hf")

# loading microsat data
data(mouse_msats)
msats <- convert_raw(mouse_msats)

test_that("g2 results matrix is computed", {
    expect_equal(is.matrix(simulate_r2_hf(n_ind = 20, subsets = c(2,4), reps = 10, mean_MLH = 0.5, 
                                       sd_MLH = 0.025)$estMat), TRUE)
    expect_equal(is.matrix(simulate_r2_hf(n_ind = 20, subsets = c(2,4), reps = 10, 
                                       genotypes = msats, mean_MLH = 0.5, 
                                       sd_MLH = 0.025)$estMat), TRUE)
    expect_equal(is.matrix(simulate_r2_hf(n_ind = 20, subsets = c(2,4), reps = 10, 
                                       genotypes = msats)$estMat), TRUE)
})

test_that("g2 results matrix is computed with snps", {
    expect_equal(is.matrix(simulate_r2_hf(n_ind = 20, subsets = c(2,4), reps = 10, mean_MLH = 0.5, 
                                       sd_MLH = 0.025, type = "snps")$estMat), TRUE)
    expect_equal(is.matrix(simulate_r2_hf(n_ind = 20, subsets = c(2,4), type = "snps", reps = 10, 
                                       genotypes = msats, mean_MLH = 0.5, 
                                       sd_MLH = 0.025)$estMat), TRUE)
    expect_equal(is.matrix(simulate_r2_hf(n_ind = 20, subsets = c(2,4), type = "snps", reps = 10, 
                                       genotypes = msats)$estMat), TRUE)
})

test_that("subsets have the right range", {
    expect_error(simulate_r2_hf(n_ind = 20, subsets = 0, reps = 10, mean_MLH = 0.5, 
                             sd_MLH = 0.025))
    expect_error(simulate_r2_hf(n_ind = 20, subsets = c(1,2), reps = 10, mean_MLH = 0.5, 
                             sd_MLH = 0.025))
    expect_error(simulate_r2_hf(n_ind = 20, subsets = NULL, reps = 10, mean_MLH = 0.5, 
                             sd_MLH = 0.025))
})

test_that("MLH have the right range", {
    expect_error(simulate_r2_hf(n_ind = 20, subsets = 0, reps = 10, mean_MLH = 1.5, 
                             sd_MLH = 0.025))
    expect_error(simulate_r2_hf(n_ind = 20, subsets = c(1,2), reps = 10, mean_MLH = -0.5, 
                             sd_MLH = 0.025))
})