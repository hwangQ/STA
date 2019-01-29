library(lpSolveAPI)
library(SuppDists)
library(catR)

##########################################################################
# Read source code
src.dir <- "E:/Simulation/6.ShadowCAT/Rcode/R/"
src.files <- list.files(src.dir, pattern="*.r")
for(files in src.files) source(paste0(src.dir, files))

# set directory
mydir <- "E:/2018 summer intership/Analysis/Rcode/Temp"
setwd(mydir)

# read item bank and content category information
bank_info <- read.csv("itempar_bank_2PL.csv")
names(bank_info) <- c("ID", "PARAM.1", "PARAM.2", "PARAM.3", "CLASS")
content <- read.csv("Content_Requirement.csv")
colnames(content) <- c("CLASS", "MIN", "MAX")
content$MAX <- c(20, 12, 8, 6)

# transform a content variable to numeric variable
x <- bank_info[, 5]
x <- as.numeric(x)
bank_info$CLASS_INT <- x

# assign other information (e.g., model, score categories) to each item in the bank
bank_info$MODEL <- "3PLM"
bank_info$CATEGORY <- 2

# reorder columns
bank_info <- bank_info[, c(1, 7, 8, 5, 6, 2, 3, 4)]

# set condtions
item.pool <- bank_info[bank_info$CATEGORY == 2, ]
lp.control <- list(timeout=5, epsint=0.1, mip.gap=c(0.1, 0.05))
true.theta <- 0
ini.theta <- 0
min.test <- 20
max.test <- 40
b.content <- list(lowb.content=content$MIN, upb.content=content$MAX)
D <- 1.702
a.strat <- TRUE
n.strat <- 3
breaks <- c(8, 8, 4)
ransq <- TRUE
par.ransq <- c(3, 2, 1)
post.ransq <- 2
inter <- list(method="WL", priorDist="norm", priorPar=c(0, 1), parInt = c(-4, 4, 81), range=c(-4, 4))
final <- list(method="ML", priorDist="norm", priorPar=c(0, 1), parInt = c(-4, 4, 81), range=c(-4, 4))
estimation <- list(inter=inter, final=final)
crit.se <- 0.2
true.thetas <- rnorm(500, 0, 1)
ini.thetas <- 0
parallel <- TRUE
useSE <- FALSE
sd.con <- 0.5
seed <- 12

##############################################################################################

ShadowCAT(item.pool=item.pool, crit.se=crit.se, min.test=min.test, max.test=max.test, true.theta=true.theta, 
		ini.theta=ini.theta, b.content=b.content, parallel=parallel, useSE=useSE, sd.con=sd.con, 
					a.strat=a.strat, n.strat=n.strat, breaks=NULL, ransq=TRUE, par.ransq, post.ransq, 
					estimation=list(inter, final), D, lp.control) 


# a function to implement a Shadow test with variable-length for a person

ShadowCAT <- function(item.pool, crit.se, min.test, max.test, true.theta, ini.theta=0, 
					b.content=list(lowb.content, upb.content), parallel=TRUE, useSE=FALSE, sd.con=0.5, 
					a.strat=TRUE, n.strat, breaks=NULL, ransq=TRUE, par.ransq, post.ransq, 
					estimation=list(inter, final), D, lp.control) {

	# decompose arguments
	lowb.content <- b.content$lowb.content
	upb.content <- b.content$upb.content
	inter.method <- estimation$inter
	final.method <- estimation$final
	
	
	# implement a minimum length of Shadow CAT
	if(parallel) {
		x <- ShadowPreSim_Parallel(item.pool=item.pool, min.test=min.test, true.theta=true.theta, ini.theta=ini.theta, 
								lowb.content=lowb.content, useSE=useSE, sd.con=sd.con, a.strat=a.strat, n.strat=n.strat, 
								breaks=breaks, ransq=ransq, par.ransq=par.ransq, inter.method=inter.method, D=D, 
								lp.control=lp.control) 
	} else {
		x <- ShadowPreSim(item.pool=item.pool, min.test=min.test, true.theta=true.theta, ini.theta=ini.theta, 
						lowb.content=lowb.content, a.strat=a.strat, n.strat=n.strat, breaks=breaks,
						ransq=ransq, par.ransq=par.ransq, inter.method=inter.method, D=D, lp.control=lp.control) 
	}
	
	# implement a post Shadow CAT with variable-length
	y <- ShadowPostSim(x=x, item.pool=item.pool, crit.se=crit.se, max.test=max.test, upb.content=upb.content, 
					post.ransq=post.ransq, inter.method=inter.method, final.method=final.method, D=D, 
					lp.control=lp.control)


	rr <- list(preShadow=x, postShadow=y)
	class(rr) <- "ShadowCAT"
	
	rr

}


test <- ShadowCATall(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas, ini.thetas=ini.thetas, 
			min.test=min.test, max.test=max.test, b.content=b.content, parallel=parallel, useSE=useSE, sd.con=sd.con,
			a.strat=a.strat, n.strat=n.strat, breaks=breaks, ransq=ransq, par.ransq=par.ransq, post.ransq=post.ransq, 
			estimation=estimation, D=D, lp.control=lp.control, seed=seed)

# a function to implement a Shadow test with variable-length for all examinees

ShadowCATall <- function(item.pool, crit.se, true.thetas, ini.thetas, min.test, max.test, b.content, 
						parallel=TRUE, useSE=FALSE, sd.con=0.5, a.strat=TRUE, n.strat, breaks=NULL, 
						ransq=TRUE, par.ransq, post.ransq, estimation, D, lp.control, seed=123) {

	# set seed number
	set.seed(seed)

	# check the number of examinees
	n.std <- length(true.thetas)
	
	# check the initial thetas
	if(length(ini.thetas) == 1L) ini.thetas <- rep(ini.thetas, n.std)
	
	# initialize the progress bar
	pb <- winProgressBar(title="Shadow CAT Progress: 0%", label="No exxaminee is done", min=0, max=100, initial=0)
	
	# make an empty list to contain the Shadow CAT results for all examinees
	catR <- vector('list', n.std)
	names(catR) <- paste0("examinee.", 1:n.std)
	for(r in 1:n.std) {

		catR[[r]] <- ShadowCAT(item.pool=item.pool, crit.se=crit.se, min.test=min.test, max.test=max.test, 
									true.theta=true.thetas[r], ini.theta=ini.thetas[r], b.content=b.content, 
									parallel=parallel, useSE=useSE, sd.con=sd.con, a.strat=a.strat, n.strat=n.strat, 
									breaks=breaks, ransq=ransq, par.ransq=par.ransq, post.ransq=post.ransq, 
									estimation=estimation, D=D, lp.control=lp.control)
		catR[[r]] <- unclass(catR[[r]])
		
		# modify the progress bar (using setWinProgressBar function)
		info <- sprintf("Shadow CAT Progress: %d%%", round((r/n.std)*100))
		setWinProgressBar(pb, r/n.std*100, title=info, label=paste0("Examinee ", r, " is done"))

	}
	
	# make a summary for final results
	n.itemT <- sapply(1:n.std, function(i) catR[[i]]$postShadow$final$n.item)
	true.theatT <- sapply(1:n.std, function(i) catR[[i]]$postShadow$metaInfo$true.theta)
	theta.estT <- sapply(1:n.std, function(i) catR[[i]]$postShadow$final$theta.est)
	se.estT <- sapply(1:n.std, function(i) catR[[i]]$postShadow$final$se.est)
	final_df <- data.frame(Length=n.itemT, ThetaTrue=true.theatT, ThetaEst=theta.estT, SeEst=se.estT)
	rownames(final_df) <- paste0("std.", 1:n.std)
	
	# make a summary for all records
	theta_recod <- lapply(1:n.std, function(i) catR[[i]]$postShadow$record$theta)
	theta_recod <- t(do.call(rowr::'cbind.fill', c(theta_recod, fill=NA)))
	se_recod <- lapply(1:n.std, function(i) catR[[i]]$postShadow$record$se)
	se_recod <- t(do.call(rowr::'cbind.fill', c(se_recod, fill=NA)))
	item_recod <- lapply(1:n.std, function(i) catR[[i]]$postShadow$record$item)
	item_recod <- t(do.call(rowr::'cbind.fill', c(item_recod, fill=NA)))
	item_ans <- lapply(1:n.std, function(i) catR[[i]]$postShadow$record$response)
	item_ans <- t(do.call(rowr::'cbind.fill', c(item_ans, fill=NA)))
	rownames(theta_recod) <- paste0("std.", 1:n.std)
	rownames(se_recod) <- paste0("std.", 1:n.std)
	rownames(item_recod) <- paste0("std.", 1:n.std)
	rownames(item_ans) <- paste0("std.", 1:n.std)
	
	summary <- list(final_df=final_df, record=list(n.std=n.std, theta=theta_recod, se=se_recod, item=item_recod, response=item_ans))
	
	###### Closing the Progress Bar
	close(pb)
	
	# return the results
	rr <- list(catR=catR, summary=summary)
	class(rr) <- "ShadowCATall"

	rr

}

# a function to implement a Shadow test with variable-length for all examinees
true.vec <- seq(-3.0, 3.0, 0.5)
nrep <- 500


for(k in 1:length(true.vec)) {
	assign(paste0("cond_true_", k), NA)
	for(r in 1:nrep) { 
		assign(paste0("cond_true_", k), 
				ShadowCATall(item.pool=item.pool, crit.se=crit.se, true.thetas=true.vec[k], ini.thetas=ini.thetas, 
							min.test=min.test, max.test=max.test, b.content=b.content, parallel=parallel, useSE=useSE, sd.con=sd.con,
							a.strat=a.strat, n.strat=n.strat, breaks=breaks, ransq=ransq, par.ransq=par.ransq, post.ransq=post.ransq, 
							estimation=estimation, D=D, lp.control=lp.control, seed=seed)
				)
	}
}



test <- ShadowCATall(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas, ini.thetas=ini.thetas, 
			min.test=min.test, max.test=max.test, b.content=b.content, parallel=parallel, useSE=useSE, sd.con=sd.con,
			a.strat=a.strat, n.strat=n.strat, breaks=breaks, ransq=ransq, par.ransq=par.ransq, post.ransq=post.ransq, 
			estimation=estimation, D=D, lp.control=lp.control, seed=seed)



## Evaluation

sapply(1:n.std, function(i) x$catR[[i]]$preShadow$record$item)

x$catR[[i]]$postShadow$final$n.item

exposRate <- function(x, item.poop) {

	nitem.pool <- nrow(item.pool)
	item_recod <- x$summary$record$item
	item_recod_f <- factor(as.numeric(item_recod), levels=1:nitem.pool)
	n.std <- x$summary$record$n.std
	er <- table(item_recod_f) / n.std
	er

}


cat_eff <- function(x) {

	n.itemT <- x$summary$final_df$Length
	n.std <- x$summary$record$n.std
	max.test <- x$catR[[1]]$postShadow$metaInfo$max.test

	mean.length <- mean(n.itemT)
	n.unreached <- sum(n.itemT == max.test) / n.std
	rr <- c(mean.length = mean.length, n.unreached=n.unreached)
	rr
	
}

BiasRmse <- function(x) {

	theta_true <- x$summary$final_df$ThetaTrue
	theta_est <- x$summary$final_df$ThetaEst
	n.std <- x$summary$record$n.std

	bias <- sum(theta_true - theta_est) / n.std
	rmse <- sqrt(sum((theta_true - theta_est)^2) / n.std)
	rr <- c(bias=bias, rmse=rmse)
	rr

}

	n.itemT <- x$summary$final_df$Length
	max.test <- x$catR[[i]]$postShadow$metaInfo$max.test

	mean.length <- mean(n.itemT)
	n.unreached <- sum(n.itemT == max.test) / n.std
	rr <- c(mean.length = mean.length, n.unreached=n.unreached)
	rr



	theta_true <- x$summary$final_df$ThetaTrue
	theta_est <- x$summary$final_df$ThetaEst
	n.std <- 3000

	bias <- sum(theta_true - theta_est) / n.std
	rmse <- sqrt(sum((theta_true - theta_est)^2) / n.std)
	rr <- c(bias=bias, rmse=rmse)
	rr