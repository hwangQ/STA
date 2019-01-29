# A function for generating response data for a test 
simdat <- function(theta, a.dc, b.dc, g.dc, a.py, b.py, cats, pModel, D) { 
	
	# check whether argument is correctly specified
	if(missing(cats)) stop("Category of each item is missing") 
	
	# Set conditions
	nstd <- length(theta)
	nitem <- length(cats)
	
	# Create an empty matrix
	simdat <- matrix(999, nrow=nstd, ncol=nitem)

	# Set initial numbers
	count.dc <- 0
	count.py <- 0

	# Simulate data
	for(i in 1:nitem) {
		
		# for dichotomous items
		if(cats[i] == 2) {
			count.dc <- count.dc + 1
			simdat[, i] <- simdat_dc(theta, a=a.dc[count.dc], b=b.dc[count.dc], g.dc[count.dc], D=D)
		}
		
		# for polytomous items
		if(cats[i] > 2) {
			count.py <- count.py + 1
			simdat[, i] <- simdat_py(theta, a=a.py[count.py], d=b.py[[count.py]], D=D, pModel=pModel[count.py])
		}
	
	}
	
	return(simdat)

}

# A function for generating binary data for one item 
simdat_dc <- function(theta, a, b, g, D) {

	# Number of examinees
	nstd <- length(theta)

	# Calculate true probability for each category
	fitted <- pl.fn(theta, a, b, g, D)

	# Sample random variables from uniform dist
	rv_unif <- runif(nstd, 0, 1)

	# Simulated Response data for one item
	simdat <- ifelse(fitted >= rv_unif, 1, 0)
	
	return(simdat)
	
}

# A function for generating categorical data for one item 
simdat_py <- function(theta, a, d, D, pModel) {

	# Number of examinees
	nstd <- length(theta)

	# add zero values for the first category for GPCM
	if(pModel == "GPCM") d <- c(0, d)
	
	# Calculate true probability for each category
	fitted <- t(sapply(1:nstd, function(i) poly.fn(theta[i], a, d, D, pModel)))

	# Sample random variables from uniform dist
	rv_unif <- runif(nstd, 0, 1)
	
	# Simulated Response data for one item
	ncat <- length(d)
	cumprob <- t(apply(fitted, 1, cumsum))
	matTF <- rv_unif >= cumprob 
	simdat <- apply(matTF, 1, sum)
	
	return(simdat)
	
}