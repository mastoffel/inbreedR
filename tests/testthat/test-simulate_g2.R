context("simulate_g2")

# loading microsat data
data(mouse_msats)
msats <- convert_raw(mouse_msats)

test_that("g2 results matrix is computed", {
    expect_equal(is.matrix(simulate_g2(n_ind = 20, subsets = c(2,4), reps = 10, mean_loc_MLH = 0.5, 
                             sd_loc_MLH = 0.025)$estMat), TRUE)
    expect_equal(is.matrix(simulate_g2(n_ind = 20, subsets = c(2,4), reps = 10, 
                                       genotypes = msats, mean_loc_MLH = 0.5, 
                                       sd_loc_MLH = 0.025)$estMat), TRUE)
    expect_equal(is.matrix(simulate_g2(n_ind = 20, subsets = c(2,4), reps = 10, 
                                       genotypes = msats)$estMat), TRUE)
})

