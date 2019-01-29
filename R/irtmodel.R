# Functions for IRT models 

# IRT 1, 2, 3PL models
pl.fn <- function(theta, a, b, g, D) {
	Da <- D * a
	z <- Da * (theta - b)
	P <- g + (1 - g) / (1 + exp(-z))
	return(P)
}

# IRT GPC model
gpcm.fn <- function (theta, score, a, steps, D) {
	Da <- D * a
	P <- exp(sum(Da * (theta - steps[1:(score+1)])))/ 
			sum(exp(cumsum(Da * (theta - steps))))
	return(P)
}

# IRT GRM model
grm.fn <- function(theta, score, a, steps, D) {
	m <- length(steps)
	allP <- sapply(1:m, function(i) pl.fn(theta, a, b=steps[i], g=0, D))
	allScores <- seq(0, m)
	if(score == min(allScores)) {
		P <- 1- allP[1]
	}
	if(score > min(allScores) & score < max(allScores)) {
		P <- allP[score] - allP[score+1] 
	}
	if(score == max(allScores)) {
		P <- allP[score]
	}
	return(P)
}

# Trace Line for polytomous IRT model
poly.fn <- function(theta, a, d, D, pModel) {
	if(pModel == "GPCM") {
		P.vec <- sapply(1:length(d), function(i) gpcm.fn(theta, score=(i-1), a, steps=d, D))
	}
	if(pModel == "GRM") {
		P.vec <- sapply(1:(length(d)+1), function(i) grm.fn(theta, score=(i-1), a, steps=d, D))
	}

	return(P.vec)
}

# TCC for polytomous IRT models
poly.exp <- function(theta, a, d, D, pModel) {
	if(pModel == "GPCM") {
		sum.w <- sum(sapply(1:length(d), function(i) gpcm.fn(theta, score=(i-1), a, steps=d, D) * (i-1)))
	}
	if(pModel == "GRM") {
		sum.w <- sum(sapply(1:(length(d)+1), function(i) grm.fn(theta, score=(i-1), a, steps=d, D) * (i-1)))
	}
	return(sum.w)
}
