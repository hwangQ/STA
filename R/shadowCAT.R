# a function to implement a Shadow test with variable-length for all examinees using a parallel package
ShadowCATall_parel <- function(item.pool, crit.se, true.thetas, ini.thetas, 
                               min.test, max.test, b.content, parallel=TRUE, useSE=FALSE, sd.con,
                               a.strat=TRUE, n.strat, breaks=NULL, ransq=TRUE, par.ransq, post.ransq, 
                               estimation, D, lp.control, source.files, seed=123) {
  
  
  # set the number of cpu cores to n - 1 cores.
  numCores <- parallel::detectCores() - 1
  
  # create a parallel processesing cluster
  cl = parallel::makeCluster(numCores, type="PSOCK")
  
  # load some specific variable names into processing cluster
  parallel::clusterExport(cl, c("item.pool", "crit.se", "true.thetas", "ini.thetas", "min.test", "max.test",
                                "b.content", "parallel", "useSE", "sd.con", "a.strat", "n.strat", "breaks",
                                "par.ransq", "post.ransq", "estimation", "D", "lp.control", "seed"), envir = environment())
  
  # load source files used for the replications into processing cluster
  lapply(1:length(source.files), function(i) parallel::clusterCall(cl, fun=source, source.files[i]))
  
  # run CAT for all replications
  f <- function(i) ShadowCATall(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas[[i]], ini.thetas=ini.thetas, 
                                min.test=min.test, max.test=max.test, b.content=b.content, parallel=parallel, useSE=useSE, sd.con=sd.con,
                                a.strat=a.strat, n.strat=n.strat, breaks=breaks, ransq=ransq, par.ransq=par.ransq, post.ransq=post.ransq, 
                                estimation=estimation, D=D, lp.control=lp.control, seed=seed)
  admin <- pbapply::pblapply(X=1:length(true.thetas), FUN=f, cl=cl) # to see the progress bar
  
  # finish
  parallel::stopCluster(cl)
  
  # return results
  admin
  
}

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
	
		x <- tryCatch({ShadowPreSim_Parallel(item.pool=item.pool, min.test=min.test, true.theta=true.theta, ini.theta=ini.theta, 
						lowb.content=lowb.content, useSE=useSE, sd.con=sd.con, a.strat=a.strat, n.strat=n.strat, breaks=breaks,
						ransq=ransq, par.ransq=par.ransq, inter.method=inter.method, D=D, lp.control=lp.control)}, 
						error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
		if(is.null(x)) {
			repeat{
				x <- tryCatch({ShadowPreSim_Parallel(item.pool=item.pool, min.test=min.test, true.theta=true.theta, ini.theta=ini.theta, 
								lowb.content=lowb.content, useSE=useSE, sd.con=sd.con, a.strat=a.strat, n.strat=n.strat, breaks=breaks,
								ransq=ransq, par.ransq=par.ransq, inter.method=inter.method, D=D, lp.control=lp.control)}, 
								error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
				if(!is.null(x)) break
			}
		}
		
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



