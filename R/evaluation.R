# a function to evalute item utilization and item exposure rate
exposRate <- function(x, item.pool, crit.over=0.3, crit.less=0.02) {

	nitem.pool <- nrow(item.pool)
	n.itemT <- x$summary$final_df$Length
	item_recod <- x$summary$record$item
	item_recod_f <- factor(as.numeric(item_recod), levels=1:nitem.pool)
	n.std <- x$summary$record$n.std
	
	# item exposure rate
	er <- table(item_recod_f) / n.std
	
	# number and rate of over-exposure items 
	over_er_num <- sum(er > crit.over)
	over_er_rate <- mean(er > crit.over)

	# number and rate of less-expousre items 
	less_er_num <- sum(er < crit.less)
	less_er_rate <- mean(er < crit.less)
	
	# number and rate of unused items 
	unused_num <- sum(er == 0)
	unused_rate <- mean(er == 0)
	
	# index of item exposure rate in variable-length test
	idx.expos <- sum((er - sum(er/nitem.pool))^2) / sum(er/nitem.pool)
	
	# test overlap rate
	comItems <- function(item_recod, comb) {
		sum(item_recod[comb[1], ] %in% item_recod[comb[2], ])
	}
	overlapItems <- combn(x=c(1:n.std), m=2, FUN=comItems, item_recod=item_recod)
	denom <- sum(overlapItems) / length(overlapItems)
	numer <- sum(n.itemT) / n.std
	overlaps <- denom / numer
	
	rr <- list(expos.rate=er, over.expos=data.frame(counts=over_er_num, rate=over_er_rate), 
			less.expos=data.frame(counts=less_er_num, rate=less_er_rate), unused=data.frame(counts=unused_num, rate=unused_rate), 
			index.expos=idx.expos, overlap.rate=overlaps)
	rr

}

# a function to evalute efficiency of CAT
cat_eff <- function(x) {

	n.itemT <- x$summary$final_df$Length
	n.std <- x$summary$record$n.std
	max.test <- x$catR[[1]]$postShadow$metaInfo$max.test
	item_recod <- x$summary$record$item
	D <- x$catR$examinee.1$preShadow$metaInfo$D

	# mean test length and unreached item rates
	mean.length <- mean(n.itemT)
	unreached.rate <- sum(n.itemT == max.test) / n.std
	
	# index of efficiency
	infos <- c()
	for(i in 1:n.std) {
		items <- item_recod[i, ][!is.na(item_recod[i, ])]
		theta <- x$summary$final_df$ThetaEst[i]
		infos[i] <- test.info(x=item.pool[items, ], theta=theta, D=D)$testInfo
	}
	eff <- sum(infos) / sum(n.itemT)

	rr <- c(mean.length = mean.length, unreached.rate=unreached.rate, efficiency=eff)
	rr
	
}

# a function to evalute mesurement accuracy
cat_acc <- function(x) {

	theta_true <- x$summary$final_df$ThetaTrue
	theta_est <- x$summary$final_df$ThetaEst
	n.std <- x$summary$record$n.std

	bias <- sum(theta_true - theta_est) / n.std
	rmse <- sqrt(sum((theta_true - theta_est)^2) / n.std)
	rr <- c(bias=bias, rmse=rmse)
	rr

}