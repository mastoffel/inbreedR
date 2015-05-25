context("convert_raw")

library(inbreedR)
data <- data("seal_microsats")

test_that("gentoypes are converted correctly", {
    expect_equal(ncol(convert_raw(seal_microsats, miss = NA)), ncol(seal_microsats)/2)
})