# a function to creat a plot for test information functions for all modules

plot.ata_mst <- function(x, xlab="Theta", ylab="TIF", main="Simultaneous Assembley", col=1:4, lwd=2, ...) {
	
	# decomposition of an object
	nstg <- x$metainfo$n.stage
	nmod <- x$metainfo$n.module
	target_rt_sc <- x$target.sc$router
	target_stg2_sc <- x$target.sc$stage2
	info_rt <- x$info.obs$router
	info_stg2 <- x$info.obs$stage2
	theta <- x$theta
	
	# target information
	ylim.max <- max(c(unlist(x$target.sc), unlist(x$info.obs))) + .5
	plot(target_rt_sc ~ theta, type='l', col=col[1], lty=2, lwd=lwd, 
		ylim=c(0, ylim.max),
		ylab=ylab, xlab=xlab, 
		main=main, ...)
	for(i in 1:nmod[1]) {
		lines(target_stg2_sc[[i]] ~ theta, lty=2, lwd=lwd, col=col[1+i])
	}

	# TIF based on selected items
	lines(info_rt ~ theta, col=col[1], lwd=lwd)
	for(i in 1:nmod[1]) {
		lines(info_stg2[[i]] ~ theta, col=col[i+1], lwd=lwd)
	}
	
	# make a legend
	legend("topleft", legend=c("Router", paste0("Stage2 (M", 1:nmod[1], ")")), title="Module", lty=1, col=col)
	legend("topright", legend=c("Router", paste0("Stage2 (M", 1:nmod[1], ")")), title="Target", lty=2, col=col)

}