
# g2_snps <- function(genotypes) {

genotypes2 <- read.table("raw_41loci_ordered.txt", row.names = 1)
genotypes <- read.table("mice_snps.txt")

g2_snps <- function(genotypes) {
# turn data into 0 (homozygote), 1 (heterozygote) or -1 (NA on one locus)
checkhet <- function(x) {

        s1 <- seq(1, length(x), 2)
        newx <- as.vector(rep(NA, length(x)/2))
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
origin <- apply(genotypes, 1, checkhet)

# define matrix with 1 for missing and 0 for all others
m <- origin
m[m==1] <- 0
m[m==-1] <- 1

# H matrix with 0 for -1
h <- origin
h[h==-1] <- 0

n <- ncol(origin) # number of individuals
l <- nrow(origin) # number of loci

# mij: proportion of individuals missing on i and j ?s locus
mtemp <- m %*% t(m)
m_ij <- mtemp/n

# vector with rowsums for missing data matrix
m_loc_temp <- apply(m, 1, sum)
m_loc <- m_loc_temp / n

# rowsum h
rowsum_h <- apply(h, 1, sum)

# variance of colsums
B_hat <- var((apply(h, 2, sum)))

#
sum_rowsum_squared <- sum(apply(h, 1, sum)^2)
C_hat <- (1/(n-1)) * (sum(h) - (1/n) * sum_rowsum_squared)

#
sum_colsum_squared <- sum(apply(h, 2, sum)^2)
A_hat <- (1/(n * (n-1))) * (sum(h)^2 - sum_rowsum_squared -
           sum_colsum_squared - sum(h))

#
M_ind <- vector()
for (i in 1:n) {
        M_temp <- vector()
        for (k in 1:l) {
                M_temp[k] <- (rowsum_h[k] * m[k, i]) / (1 - m_loc[k])
        }
        M_ind[i] <- (sum(M_temp))^2 - sum(M_temp^2)
}

M_hat <- (1/(n^2)) * sum(M_ind)

#
X_temp <- vector()
for (k in 1:l) {
        X_temp[k] <- (rowsum_h[k] * m_loc[k]) / (1 - m_loc[k])
}

a_hat_temp <- (1/n) * (sum(X_temp^2) - (sum(X_temp))^2)

a_hat_nenn <- (1/n) * ((sum(h))^2 - sum_rowsum_squared)

a_hat <- (M_hat + a_hat_temp) / a_hat_nenn

g2_snps <- ((1 + (B_hat - C_hat)/A_hat) / (1 + a_hat)) - 1

}









