#' Simulate g2 
#' 
#' This function can be used to simulate genotype data, draw random subsamples and calculate the
#' respective g2 values. Every subsample of markers is drawn independently to give insights
#' into the variation of g2 values calculated from a given number of markers and individuals. 
#' Empirical genotype data can be used to estimate multilocus heterozygosity (MLH) and provide a distribution
#' (from mean and sd of MLH) for the simulated genotypes.
#'
#' @param n_ind number of individuals to sample from the population
#' @param subsets a vector specifying the sizes of marker-subsets to draw. For a subset of 20 markers, subsets = c(2, 5, 10, 15, 20) could
#'        be a reasonable choice. The minimum subset size is 2 and the maximum is the number of markers in the data
#' @param reps number of resampling repetitions
#' @param type specifies g2 formula to take. Type "snps" for large datasets and "msats" for smaller datasets.
#' @param genotypes optional: provide genotypes in \code{inbreedR} format to estimate the empirical multilocus heterozygosity
#'        (MLH) and use these estimates to simulate data from the same distribution
#' @param mean_loc_MLH mean multilocus heterozygosity. Should range between 0-1. This will be automatically
#'        calculated from the empirical data when genotypes in \code{inbreedR} format are given.
#' @param sd_loc_MLH standard deviation of multilocus heterozygosity. This will be automatically
#'        calculated from the empirical data when genotypes in \code{inbreedR} format are given.
#' @param CI Confidence intervals to calculate (default to 0.95)
#'
#' @details The \code{simulate_g2} function can be used to explore the confidence of g2 estimates calculated
#'          from marker sets of different sizes. Every new locus set is drawn independently.
#'          The simulation assumes (1) unlinked loci, (2) equal allele frequencies among all loci 
#'          with expected heterozygosity of 0.5 and (3) random mating. The mean and standard 
#'          deviation of genome-wide expected heterozygosity can either be specified with 
#'          the \code{mean_loc_MLH} and \code{sd_loc_MLH} arguments or will be calculated 
#'          from emprical genotypes when these are specified in the \code{genotypes} argument.
#'          
#'          
#' @author  Marty Kardos (marty.kardos@@ebc.uu.se) &
#'          Martin A. Stoffel (martin.adam.stoffel@@gmail.com) 
#'          
#' @examples 
#' data(mouse_msats)
#' genotypes <- convert_raw(mouse_msats)
#' sim_g2 <- simulate_g2(n_ind = 10, subsets = c(2,4,6), reps = 100, genotypes = genotypes)
#' plot(sim_g2)
#' @export


simulate_g2 <- function(n_ind = NULL, subsets = NULL, reps = 100, type = c("msats", "snps"),
                        genotypes = NULL, mean_loc_MLH = NULL, sd_loc_MLH = NULL,
                        CI = 0.95) {
################################################################################
# simulate a population with variable inbreeding
# then estimate g2 from independently sampled / non-overlapping subsets of loci 
################################################################################

# predefine
# check if empirical data is given
if (!is.null(genotypes)) {
    if (check_data(genotypes)) {
        loc_MLH <- MLH(genotypes)
        mean_loc_MLH <- mean(loc_MLH)
        sd_loc_MLH <- sd(loc_MLH)
    }
}

# number of individuals to sample from the population
if (is.null(n_ind))   stop("Specify the number of individuals to sample with n_ind")  
# subsets of loci to sample   
if (is.null(subsets)) stop("specify the size of loci subsamples in 'subsets', i.e. subsets = c(2,4,6,8) to 
                           calculate g2 from up to 8 loci")
if (any(subsets < 2)) stop("Specify a minimum of 2 markers in subsets")
if (!isTRUE(all(subsets == floor(subsets)))) stop("'subsets' must only contain integer values")
# expected heterozygosity at all loci (assumes expected heterozygosity has zero variance across a the simulated loci)
if (is.null(mean_loc_MLH)) stop("Specify the mean expected heterozygosity (value between 0-1) with mean_loc_MLH")     
if (is.null(sd_loc_MLH))  stop("Specify the standard deviation of expected heterozygosity with sd_loc_MLH")

    # check g2 function argument
    if (length(type) == 2){
        type <- "msats"
    } else if (!((type == "msats")|(type == "snps"))){
        stop("type argument needs to be msats or snps")
    } 
    
    # define g2 function
    if (type == "msats") {
        g2_fun <- g2_microsats
    } else if (type == "snps") {
        g2_fun <- g2_snps
    }
    
# total number of loci to simulate
n_loc <- subsets[length(subsets)]
allLoci <- reps*n_loc                            

#-------------------------
# simulate the individual 
# genotypes
#-------------------------

hets <- NULL        # initialize a data frame to store the genotypic information

for (i in 1:n_ind) {
    
	thisHet <- NULL                       # randomly select a TRUE genome-wide MLH (i.e., the proportion of hypothetically infinitely many loci that are heterozygous in the ith individual)
	thisHet <- rnorm(1,mean=mean_loc_MLH,sd=sd_loc_MLH)

	rands <- NULL                         # randomly generated numbers between 0 and 1 that are used to determine whether the individual is heterozygous at each locus
	rands <- runif(allLoci,min=0,max=1)

	theseHets <- NULL
	theseHets <- as.numeric(rands < thisHet)
	
	hets <- rbind(hets,theseHets)
	}

#------------------------------------------------------------------------
# repetitively subsample the loci independently, each time estimating g2
#------------------------------------------------------------------------
sampNVec <- c(subsets)

estMat <- NULL    # matrix to store teh g2 estimates
sampCols <- 1:ncol(hets)    # vector of loci that are available for sampling


estMat <- NULL

for (i in 1:length(sampNVec))    # loop through the differen subsample sizes
	{
	theseEsts <- rep(NA,reps)    # vector to store the estimates from this number of loci
	for (j in 1:reps)
		{
        theseSampCols <- NULL                                       # get a new independent sample of loci
		theseSampCols <- sample(sampCols,sampNVec[i],replace=FALSE)

		theseGenos <- NULL
		theseGenos <- hets[,theseSampCols]

		theseEsts[j] <- g2_fun(theseGenos)[2][[1]]
		sampCols <- sampCols[-theseSampCols]

		}

	estMat <- rbind(estMat,theseEsts)
	sampCols <- 1:ncol(hets)    # reconstitute the original vector of loci that are available for sampling
	print(paste("done with subsampling of ",sampNVec[i]," loci",sep=""))

	}
estMat <- unname(estMat)
# get the upper and lower bounds of the y-axis

minG2 <- min(estMat)
maxG2 <- max(estMat)

# calculate CIs and SDs
calc_CI <- function(estMat_subset) {
    stats::quantile(estMat_subset, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
}

all_CI <- t(apply(estMat, 1, calc_CI))
all_sd <- apply(estMat, 1, sd)

res <- list(call=match.call(),
            estMat = estMat,
            n_ind = n_ind,
            n_loc = n_loc,
            subsets = subsets,
            reps = reps,
            genotypes = genotypes,
            mean_loc_MLH = mean_loc_MLH,
            sd_loc_MLH = sd_loc_MLH,
            minG2 = minG2,
            maxG2 = maxG2,
            sampNVec = sampNVec,
            all_CI = all_CI,
            all_sd = all_sd
            )

class(res) <- "inbreed"
return(res)
}


# additional  function to calculate marker MLH for simulation
MLH <- function(genotypes) {
    # transform to matrix
    genes <- as.matrix(genotypes)
    # genes[is.na(genes)] <- -1
    # get logical matrix of non-missing values as TRUE
    typed <- (genes == 1) | (genes == 0)
    typed[is.na(typed)] <- FALSE
    # initialise
    nloc <- ncol(genes)
    nind <- nrow(genes)
    typed_sum <- colSums(typed, na.rm = TRUE)
    # heterozygosity per locus
    het_loc <- colSums(genes == 1, na.rm = TRUE) / typed_sum
    # replicate vector to matrix
    het_loc_mat <- matrix(het_loc, nrow = nind, ncol = nloc, byrow = TRUE)
    het_loc_mat[!typed] <- 0
    mh <- rowSums(het_loc_mat, na.rm = TRUE)
    N  <- rowSums(typed, na.rm = TRUE)
    H  <- rowSums(genes == 1, na.rm = TRUE)
    MLH <- (H/N)
}




