# turn raw genotype data into 0 (homozygote), 1 (heterozygote) or -1 (NA on one locus)

convert_raw <- function(genotypes, NAval = NA) {

genotypes[genotypes == NAval] <- NA

check_het <- function(x) {

        if ((length(x) %% 2) == 1) {
                s1 <- seq(1, (length(x)-1), 2)
                newx <- as.vector(rep(NA, ceiling(length(x)/2)))
                }
        if ((length(x) %% 2) == 0) {
                s1 <- seq(1, length(x), 2)
                newx <- as.vector(rep(NA, length(x)/2))
                }

        count <- 1

        for(i in s1){
                if (is.na(x[i] | x[i + 1])) {
                        newx[count] = -1
                        count = count + 1
                } else if (x[i] == x[i + 1]) {
                        newx[count] = 0
                        count = count + 1
                } else if (x[i] != x[i + 1]) {
                        newx[count] = 1
                        count = count + 1
                }
        }
        newx
}

# original full data matrix
origin <- as.data.frame(t(apply(genotypes, 1, check_het)))
row.names(origin) <- row.names(genotypes)
origin
}
