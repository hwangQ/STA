# a function to assign a path to each router score
give_path <- function(score, cut.score) { 

	# number of scores
	nscore <- length(score)

	# number of categories (or levels)
	cats <- length(cut.score) + 1

	# set total lables to be given to each examinee
	level <- 1:cats

	# create path vectors
	path <- rep(level[1], nscore)

	# give a level correspoding each score
	for(i in 1:length(cut.score)) {
		path[which(score >= cut.score[i])] <- level[i+1]
	}

	rr <- list(path=path, score=score, cut.score=cut.score)
	rr

}

# calculate moments
cal_moments <- function(node, weight) {

	mu <- sum(node * weight)
	sigma2 <- sum(node^2 * weight) - mu^2
	rr <- c(mu=mu, sigma2=sigma2)
	rr

}

# marginal reliability
mrel <- function(var.pop, var.cond, wieghted=TRUE, w) {
	
	if(wieghted) {
		rr <- (var.pop - weighted.mean(x=var.cond, w=w)) / var.pop
	} else {
		rr <- (var.pop - mean(x=var.cond)) / var.pop
	}
	
	rr

}


# calculation of four moments of normal distribution
norm_moments <- function(mu, sigma) {

  m1 <- mu
  m2 <- mu^2 + sigma^2
  m3 <- mu * (mu^2 + 3 * sigma^2)
  m4 <- mu^4 + 6 * mu^2 * sigma^2 + 3 * sigma^4
  
  rr <- c(m1=m1, m2=m2, m3=m3, m4=m4)
  rr
  
}