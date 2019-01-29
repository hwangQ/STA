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
content$MAX <- c(20, 10, 6, 4)

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
true.theta <- 1
ini.theta <- 0
min.test <- 20
max.test <- 40
lowb.content <- content$MIN
upb.content <- content$MAX
D <- 1.702
a.strat <- TRUE
n.strat <- 3
breaks <- c(8, 8, 4)
ransq <- TRUE
par.ransq <- c(3, 2, 1)
post.ransq <- 2
inter.method <- list(method="WL", priorDist="norm", priorPar=c(0, 1), parInt = c(-4, 4, 81), range=c(-4, 4))
final.method <- list(method="ML", priorDist="norm", priorPar=c(0, 1), parInt = c(-4, 4, 81), range=c(-4, 4))
crit.se <- 0.2
useSE <- FALSE
sd.con <- 0.5

##################################################################
# A function to simulate a minimum length of Shadow CAT

ShadowPreSim_Parallel <- function(item.pool, min.test, true.theta, ini.theta=0, lowb.content, useSE=FALSE, sd.con=0.5,
					a.strat=TRUE, n.strat, breaks=NULL, ransq=TRUE, par.ransq, inter.method, D, lp.control) {


	# check the appropriateness of input of arguments
	if(min.test != sum(lowb.content)) stop("A length of minimum test should be the same as the sum of low bounds of content")
	
	###########################################################################
	### 1) PREPERATAION FOR SHADOW CAT
	###########################################################################
	# set conditions
	F <- ifelse(useSE, 3, 2) 
	cats <- item.pool$CATEGORY
	model <- item.pool$MODEL
	item.id <- item.pool$ID
	a <- item.pool$PARAM.1; b <- item.pool$PARAM.2; g <- item.pool$PARAM.3
	cont <- item.pool$CLASS_INT
	C <- length(unique(cont)) # number of content categories
	I <- length(cont) # total number of items in a pool
	NcMin <- lowb.content
	if(a.strat) {
		if(is.null(breaks)) {
			Na <- rep(min.test/n.strat, n.strat) 
			NaCum <- cumsum(Na)
		}
		if(length(breaks) == 1L) {
			Na <- rep(min.test/breaks, breaks) 
			NaCum <- cumsum(Na)
		}
		if(length(breaks) > 1L) {
			Na <- breaks 
			NaCum <- cumsum(Na)
		}
	}
	
	# assign item to each category of content
	Vc <- list()
	for(i in 1:C) {
		Vc[[i]] <- c(1:I)[cont == i]
	}
	
	# straticiation item pool based on a parameters
	a_rank <- rank(item.pool$PARAM.1)
	idx_astr <- as.numeric(cut(a_rank, breaks=n.strat, labels=1:n.strat))

	# assign item to each category of a-stratification
	if(a.strat) {
		Va <- list()
		for(i in 1:n.strat) {
			Va[[i]] <- c(1:I)[idx_astr == i]
		}
	}

	# set the first ability to construct a first shadow test
	theta.est <- ini.theta
	
	# create empty vectors that the CAT procedures are recorded
	item_recod <- c(); item_ans <- c(); theta_recod <- c(); se_recod <- c()
	shadow_recod <- array(NA, c(min.test, min.test))

	###########################################################################
	### 3) MAKE A SHADOW TEST
	###########################################################################
	# make a model
	shadow <- lpSolveAPI::make.lp(nrow=0, ncol=I * F)

	# set control parameters: maximization problem; integer tolerance is et to 0.1; 
	# absolute MIP gaps is set to 0.1, relative MIP gap is set to 0.5
	lpSolveAPI::lp.control(lprec=shadow, sense="max", timeout=lp.control$timeout, epsint=lp.control$epsint, mip.gap=lp.control$mip.gap)

	# contratint: type of decision variables
	lpSolveAPI::set.type(lprec=shadow, columns=1:(I*F), type="binary")
	lpSolveAPI::set.bounds(lprec=shadow, lower=rep(0, I*F), upper=rep(1, I*F))

	# constraint: test length
	for(f in 1:F) {
		lpSolveAPI::add.constraint(lprec=shadow, xt=rep(1, I), type="=", rhs=min.test, indices=(I*(f-1)+1):(I*f))
	}
	
	# constraints: content
	for(f in 1:F) {
		for(i in 1:length(NcMin)) {
			lpSolveAPI::add.constraint(lprec=shadow, xt=rep(1, length(Vc[[i]])), type=">=", rhs=NcMin[i], indices=(f*I-I) + Vc[[i]])
		}
	}
	
	if(a.strat) {
		# constraints: a-stratification
		for(f in 1:F) {
			for(i in 1:(length(Na)-1)) {
				lpSolveAPI::add.constraint(lprec=shadow, xt=rep(1, length(Va[[i]])), type=">=", rhs=Na[i], indices=(f*I-I) + Va[[i]])
			}
		}
	}
	
	for(i in 1:min.test) {

		###################################
		## (1) create a Shadow test with a current ability estimate
		###################################
		# calcuate item information for every items in a pool based on a current ability estimates
		if(useSE) {
			if(i == 1L) {
				info_1 <- as.numeric(test.info(x=item.pool, theta=theta.est, D=D)$itemInfo)
				info_2 <- as.numeric(test.info(x=item.pool, theta=(theta.est - .5), D=D)$itemInfo)
				info_3 <- as.numeric(test.info(x=item.pool, theta=(theta.est + .5), D=D)$itemInfo)
			} else {
				info_1 <- as.numeric(test.info(x=item.pool, theta=theta.est, D=D)$itemInfo)
				info_2 <- as.numeric(test.info(x=item.pool, theta=(theta.est - sd.con*se.est), D=D)$itemInfo)
				info_3 <- as.numeric(test.info(x=item.pool, theta=(theta.est + sd.con*se.est), D=D)$itemInfo)
			}
			info <- c(info_1, info_2, info_3)
		} else {
				info <- as.numeric(test.info(x=item.pool, theta=theta.est, D=D)$itemInfo)
				info <- c(info, info)
		}
		
		
		# objective function
		lpSolveAPI::set.objfn(lprec=shadow, obj=info)
	
		# constraint: include already administered items in the previous stages
		if(i > 1) {
			for(f in 1:F) {
				lpSolveAPI::add.constraint(lprec=shadow, xt=1, type="=", rhs=1, indices=(f*I-I) + item_recod[i - 1])
			}
		}
		
		# count the number of constraints
		constr_num1 <- length(lpSolveAPI::get.constr.value(shadow))
		
		# only when a parallel shadow test forms are constructed
		if(!useSE) {
			# constraint: no overlap between two shadow tests except the already administered items
			if(i == 1L) {
				for(j in 1:I) {
					lpSolveAPI::add.constraint(lprec=shadow, xt=rep(1, 2), type="<=", rhs=1, indices=c(j, (I+j)))
				}
			} 
			if(i > 1L) {
				freeitems <- c(1:I)[!(1:I %in% item_recod)]
			for(j in freeitems) {
					lpSolveAPI::add.constraint(lprec=shadow, xt=rep(1, 2), type="<=", rhs=1, indices=c(j, (I+j)))
				}
			}
		}
		
		# count the number of constraints
		constr_num2 <- length(lpSolveAPI::get.constr.value(shadow))
		
		# solve teh model
		res_flag <- lpSolveAPI::solve.lpExtPtr(shadow)
		if(res_flag > 1) stop(paste0("A status code for constructing a Shadow test at stage ", i, " is ", res_flag))
		
		# retrieve the values of the decision variables 
		sim_opt <- lpSolveAPI::get.variables(shadow)

		# delete constraints only when the parallel shadow test forms are constructed
		if(!useSE) {
			for(j in 1:(constr_num2-constr_num1))
			lpSolveAPI::delete.constraint(lprec=shadow, constraints=(constr_num1+1))
		}
		
		###################################
		## (2) select a next item
		###################################
		# select all items consisting of a Shadow test among item pool
		sim_opt_list <- vector('list', F)
		item_shadow <- vector('list', F)
		for(f in 1:F) {
			sim_opt_list[[f]] <- sim_opt[(f*I-I+1):(f*I)]
			item_shadow[[f]] <- c(1:I)[sim_opt_list[[f]] == 1]
		}
		
		## when a-stratification is used
		if(a.strat) {

			# give a level of a-stratification to each item in the Shadow test
			if(is.null(breaks)) {
				Level <- as.numeric(cut(1:min.test, breaks=n.strat, level=1:n.strat))
			}
			if(length(breaks) == 1L) {
				Level <- as.numeric(cut(1:min.test, breaks=breaks, level=1:n.strat))
			}
			if(length(breaks) > 1L) {
				Level <- as.numeric(cut(1:min.test, breaks=c(0, NaCum), level=1:n.strat))
			}
			
			# select the candidate items that could be administered for a next item
			# items that is classified to the a-stratum corresponding to the stage are only selected
			# for example, if the corresponding a-stratum for the ith is the first stratum, then only items
			# are classified to the first a-stratum are selected among items in the Shadow test
			item_adm <- vector('list', F)
			if(Level[i] < n.strat) {
				for(f in 1:F) {
					item_adm[[f]] <- item_shadow[[f]][idx_astr[item_shadow[[f]]] == Level[i]]
				}
			} else {
				item_adm <- item_shadow
			}
	
			# exclude the items already administered in the prevous stage among the possible candidate items
			if(i > 1) {
				for(f in 1:F) {
					item_adm[[f]] <- item_adm[[f]][!item_adm[[f]] %in% item_recod]
				}
			}
			
			# item information for the seleced candidate items
			info_adm <- vector('list', F)
			for(f in 1:F) {
				info_adm[[f]] <- info[(f*I-I) + item_adm[[f]]]
			}
			
			## randomesque method is used
			item_ran <- vector('list', F)
			if(ransq) {
				
				for(f in 1:F) {
					# select the final candidate items that have high information values at the current ability estimate using a randomesque method
					item_ran[[f]] <- item_adm[[f]][rank(info_adm[[f]]) %in% 1:par.ransq[Level[i]]]
				
				}
				# randomly selet an item to be administered at the next stage
				item_next <- ifelse(length(unique(unlist(item_ran))) > 1, sample(unique(unlist(item_ran)), 1), unique(unlist(item_ran)))
				
			## randomesque method is not used
			} else {
				
				for(f in 1:F) {
					# a next item is selected only based on the maximum value of item information
					item_ran[[f]] <- item_adm[[f]][info_adm[[f]] == max(info_adm[[f]])]
				}
				# randomly selet an item to be administered at the next stage
				item_next <- ifelse(length(unique(unlist(item_ran))) > 1, sample(unique(unlist(item_ran)), 1), unique(unlist(item_ran)))
			}

		## when a-stratification is not used
		} else {

			# exclude the items already administered in the prevous stage among the possible candidate items
			item_adm <- item_shadow
			if(i > 1) {
				for(f in 1:F) {
					item_adm[[f]] <- item_adm[[f]][!item_adm[[f]] %in% item_recod]
				}
			}
			
			# a next item is selected only based on the maximum value of item information
			info_adm <- vector('list', F)
			for(f in 1:F) {
				info_adm[[f]] <- info[(f*I-I) + item_adm[[f]]]
			}
			item_ran <- vector('list', F)
			for(f in 1:F) {
				# a next item is selected only based on the maximum value of item information
				item_ran[[f]] <- item_adm[[f]][info_adm[[f]] == max(info_adm[[f]])]
			}
			item_next <- ifelse(length(unique(unlist(item_ran))) > 1, sample(unique(unlist(item_ran)), 1), unique(unlist(item_ran)))
		
		}

		# record the finally selected item
		item_recod <- c(item_recod, item_next)
	
		# recod the shadow test
		shadow_recod[i, 1:length(item_recod)] <- c(item_recod)

		# random generation of response pattern the selected item with the current ability estimate
		answer <- simdat(theta=true.theta, a.dc=a[item_next], b.dc=b[item_next], g.dc=g[item_next], cats=2, D=D)
	
		# record the generated response
		item_ans <- c(item_ans, answer)
	
		# prepare item parameters for all administered items to estimate a recent ability
		item_par <- data.matrix(item.pool[item_recod, c("PARAM.1", "PARAM.2", "PARAM.3")])
		it <- cbind(item_par, rep(1, length(item_ans)))
	
		# interim ability estimate and its standard error
		theta.est <- catR::thetaEst(it=it, x=item_ans, D=D, method=inter.method$method, priorDist=inter.method$priorDist,
							priorPar=inter.method$priorPar, parInt=inter.method$parInt, range=inter.method$range)
		se.est <- catR::semTheta(thEst=theta.est, it=it, x=item_ans, D=D, method="ML", range=inter.method$range)
	
		# recode the ability estimate and its standard error
		theta_recod <- c(theta_recod, theta.est)
		se_recod <- c(se_recod, se.est)

	}

	# give column and row names to recorded results
	col.name <- paste0("item.", 1:min.test)
	names(item_recod) <- col.name
	names(item_ans) <- col.name
	names(theta_recod) <- col.name
	names(se_recod) <- col.name
	rownames(shadow_recod) <- paste0("shadow.", 1:min.test)
	colnames(shadow_recod) <- col.name
	
	# information (e.g., item parameters, content) for all of administered items
	item_df <- item.pool[item_recod, ]
	
	#
	if(!a.strat) {
		n.strat <- NA; breaks <- NA
	}
	if(is.null(breaks)) breaks <- NA
	if(!ransq) par.ransq <- NA
	if(!useSE) sd.con <- NA
	
	rr <- list(record=list(theta=theta_recod, se=se_recod, item=item_recod, response=item_ans, shadow=shadow_recod), 
				final=list(theta.est=theta.est, se.est=se.est, test=item_df),
				metaInfo=list(true.theta=true.theta, ini.theta=ini.theta, D=D, min.test=min.test, lowb.content=lowb.content,
							useSE=useSE, sd.con=sd.con, a.strat=a.strat, n.strat=n.strat, breaks=breaks, ransq=ransq, 
							par.ransq=par.ransq),
				estimation=list(method=inter.method$method, priorDist=inter.method$priorDist, priorPar=inter.method$priorPar,
								parInt=inter.method$parInt, range=inter.method$range)
				)

	rr

}

list(method="WL", priorDist="norm", priorPar=c(0, 1), parInt = c(-4, 4, 81), range=c(-4, 4))

#################################################################################

x <- ShadowPreSim_Parallel(item.pool=item.pool, min.test=min.test, true.theta=true.theta, ini.theta=ini.theta, 
			lowb.content=lowb.content, useSE=useSE, sd.con=sd.con, a.strat=a.strat, n.strat=n.strat, breaks=breaks,
			ransq=ransq, par.ransq=par.ransq, inter.method=inter.method, D=D, lp.control=lp.control) 

x <- ShadowPreSim(item.pool=item.pool, min.test=min.test, true.theta=true.theta, ini.theta=ini.theta, 
			lowb.content=lowb.content, a.strat=a.strat, n.strat=n.strat, breaks=breaks,
			ransq=ransq, par.ransq=par.ransq, inter.method=inter.method, D=D, lp.control=lp.control) 

# A function to simulate a post Shadow CAT

ShadowPostSim <- function(x, item.pool, crit.se, max.test, upb.content, post.ransq, inter.method, final.method, D, lp.control) {


	###########################################################################
	### 1) PREPERATAION FOR SHADOW CAT
	###########################################################################
	# set conditions
	cats <- item.pool$CATEGORY
	model <- item.pool$MODEL
	item.id <- item.pool$ID
	a <- item.pool$PARAM.1; b <- item.pool$PARAM.2; g <- item.pool$PARAM.3
	cont <- item.pool$CLASS_INT
	C <- length(unique(cont)) # number of content categories
	I <- length(cont) # total number of items in a pool
	NcMax <- upb.content
	
	# extract needed information from an object
	se.est <- x$final$se.est
	true.theta <- x$metaInfo$true.theta
	theta.est <- x$final$theta.est
	theta_recod <- x$record$theta
	se_recod <- x$record$se
	item_recod <- x$record$item
	item_ans <- x$record$response
	shadow_recod_min <- x$record$shadow
	min.test <- x$metaInfo$min.test
	
	# assign item to each category of content
	Vc <- list()
	for(i in 1:C) {
		Vc[[i]] <- c(1:I)[cont == i]
	}

	# create empty vectors that the CAT procedures are recorded
	shadow_recod <- array(NA, c(max.test, max.test))
	shadow_recod[1:min.test, 1:min.test] <- shadow_recod_min

	if(se.est >= crit.se) {
	
		###########################################################################
		### 2) MAKE A SHADOW TEST
		###########################################################################
		# make a model
		shadow <- lpSolveAPI::make.lp(nrow=0, ncol=I)

		# set control parameters: maximization problem; integer tolerance is et to 0.1; 
		# absolute MIP gaps is set to 0.1, relative MIP gap is set to 0.5
		lpSolveAPI::lp.control(lprec=shadow, sense="max", timeout=lp.control$timeout, epsint=lp.control$epsint, mip.gap=lp.control$mip.gap)

		# contratint: type of decision variables
		lpSolveAPI::set.type(lprec=shadow, columns=1:I, type="binary")
		lpSolveAPI::set.bounds(lprec=shadow, lower=rep(0, I), upper=rep(1, I))

		# constraint: test length
		lpSolveAPI::add.constraint(lprec=shadow, xt=rep(1, I), type="=", rhs=max.test, indices=1:I)

		# constraints: content
		for(i in 1:length(NcMax)) {
			lpSolveAPI::add.constraint(lprec=shadow, xt=rep(1, length(Vc[[i]])), type="<=", rhs=NcMax[i], indices=Vc[[i]])
		}

		# constraint: include already administered items in the minimum length Shadow CAT
		for(i in 1:min.test) {
			lpSolveAPI::add.constraint(lprec=shadow, xt=1, type="=", rhs=1, indices=item_recod[i])
		}
	
	
		for(i in (min.test+1):max.test) {

			###################################
			## (1) create a Shadow test with a current ability estimate
			###################################
			# calcuate item information for every items in a pool based on a current ability estimates
			info <- as.numeric(test.info(x=item.pool, theta=theta.est, D=D)$itemInfo)
	
			# objective function
			lpSolveAPI::set.objfn(lprec=shadow, obj=info)
	
			# constraint: include already administered items in the minimum length Shadow CAT
			if(i > (min.test+1)) {
				lpSolveAPI::add.constraint(lprec=shadow, xt=1, type="=", rhs=1, indices=item_recod[i -1])
			}

			# solve teh model
			res_flag <- lpSolveAPI::solve.lpExtPtr(shadow)
			if(res_flag > 1) stop(paste0("A status code for constructing a Shadow test at stage ", i, " is ", res_flag))

			# retrieve the values of the decision variables 
			sim_opt <- lpSolveAPI::get.variables(shadow)

			###################################
			## (2) select a next item
			###################################
			# select all items consisting of a Shadow test among item pool
			item_shadow <- c(1:I)[sim_opt == 1]
		
			# exclude the items already administered in the prevous stage among the possible candidate items
			item_adm <- item_shadow
			item_adm <- item_adm[!item_adm %in% item_recod]
		
			# a next item is selected only based on the maximum value of item information
			info_adm <- info[item_adm]
		
			# randomesque method is used
			item_ran <- item_adm[rank(info_adm) %in% 1:post.ransq]

			# randomly selet an item to be administered at the next stage
			item_next <- ifelse(length(item_ran) > 1, sample(item_ran, 1), item_ran)
		
			# record the finally selected item
			item_recod <- c(item_recod, item_next)
	
			# recod the shadow test
			remained <- item_shadow[!item_shadow %in% item_recod]
			shadow_recod[i, ] <- c(item_recod, remained)

			# random generation of response pattern the selected item with the current ability estimate
			answer <- simdat(theta=true.theta, a.dc=a[item_next], b.dc=b[item_next], g.dc=g[item_next], cats=2, D=D)
	
			# record the generated response
			item_ans <- c(item_ans, answer)
	
			# prepare item parameters for all administered items to estimate a recent ability
			item_par <- data.matrix(item.pool[item_recod, c("PARAM.1", "PARAM.2", "PARAM.3")])
			it <- cbind(item_par, rep(1, length(item_ans)))
	
			# interim ability estimate and its standard error
			theta.est <- catR::thetaEst(it=it, x=item_ans, D=D, method=inter.method$method, priorDist=inter.method$priorDist,
								priorPar=inter.method$priorPar, parInt=inter.method$parInt, range=inter.method$range)
			se.est <- catR::semTheta(thEst=theta.est, it=it, x=item_ans, D=D, method="ML", range=inter.method$range)
	
			# recode the ability estimate and its standard error
			theta_recod <- c(theta_recod, theta.est)
			se_recod <- c(se_recod, se.est)

			# make a decision based on a critera value of SE
			if(se.est < crit.se) break
		
		}

	} else {
	
		# prepare item parameters for all administered items to estimate a recent ability
		item_par <- data.matrix(item.pool[item_recod, c("PARAM.1", "PARAM.2", "PARAM.3")])
		it <- cbind(item_par, rep(1, length(item_ans)))
	
	}
	
	# total number of administered items
	n.item <- nrow(item_par)
	
	# final estimation of ability
	theta.final <- catR::thetaEst(it=it, x=item_ans, D=D, method=final.method$method, priorDist=final.method$priorDist,
						priorPar=final.method$priorPar, parInt=final.method$parInt, range=final.method$range)
	se.final <- catR::semTheta(thEst=theta.final, it=it, x=item_ans, D=D, method=final.method$method, priorDist=final.method$priorDist,
						priorPar=final.method$priorPar, parInt=final.method$parInt, range=final.method$range)

	# give column and row names to recorded results
	col.name <- paste0("item.", 1:length(item_recod))
	names(item_recod) <- col.name
	names(item_ans) <- col.name
	names(theta_recod) <- col.name
	names(se_recod) <- col.name
	rownames(shadow_recod) <- paste0("shadow.", 1:max.test)
	colnames(shadow_recod) <- paste0("item.", 1:max.test)
	
	# information (e.g., item parameters, content) for all of administered items
	item_df <- item.pool[item_recod, ]


	rr <- list(record=list(theta=theta_recod, se=se_recod, item=item_recod, response=item_ans, shadow=shadow_recod), 
				final=list(theta.est=theta.final, se.est=se.final, test=item_df, n.item=n.item),
				metaInfo=list(crit.se=crit.se, true.theta=true.theta, D=D, min.test=min.test, max.test=max.test, upb.content=upb.content, post.ransq=post.ransq),
				estimation=list(inter=list(method=inter.method$method, priorDist=inter.method$priorDist, priorPar=inter.method$priorPar,
											parInt=inter.method$parInt, range=inter.method$range),
								final=list(method=final.method$method, priorDist=final.method$priorDist, priorPar=final.method$priorPar,
											parInt=final.method$parInt, range=final.method$range))
				)
	
	rr

}


y <- ShadowPostSim(x=x, item.pool=item.pool, crit.se=crit.se, max.test=max.test, upb.content=upb.content, 
					post.ransq=post.ransq, inter.method=inter.method, final.method=final.method, D=D, 
					lp.control=lp.control)


