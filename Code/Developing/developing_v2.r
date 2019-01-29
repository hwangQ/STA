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
cont_require <- read.csv("Content_Requirement.csv")

# transform a content variable to numeric variable
x <- bank_info[, 5]
x <- as.numeric(x)
bank_info$CLASS_RE <- x

# assign other information (e.g., model, score categories) to each item in the bank
bank_info$MODEL <- "3PLM"
bank_info$CATEGORY <- 2

# reorder columns
bank_info <- bank_info[, c(1, 7, 8, 5, 6, 2, 3, 4)]

# set condtions
lp.control <- list(timeout=5, epsint=0.1, mip.gap=c(0.1, 0.05))
true.theta <- 1
ini.theta <- 0
n.test <- 24
D <- 1.702
n.astr <- 3
ransq <- c(4, 3, 2)

##################################################################

	item.pool <- bank_info
	
	# set conditions
	cats <- item.pool$CATEGORY
	model <- item.pool$MODEL
	item.id <- item.pool$ID
	df_dc <- item.pool[item.pool$CATEGORY == 2, ]
	a <- df_dc$PARAM.1
	b <- df_dc$PARAM.2
	g <- df_dc$PARAM.3
	prm_dc <- list(a, b, g)
	cont <- item.pool$CLASS_RE

	C <- length(unique(cont)) # number of content categories
	I <- length(cont) # total number of items in a pool
	Nc <- cont_require$Router
	Na <- rep(n.test/n.astr, n.astr) 
	NaCum <- cumsum(Na)
	
	
	# assign item to each category of content
	Vc <- list()
	for(i in 1:C) {
		Vc[[i]] <- c(1:I)[cont == i]
	}
	
	# straticiation item pool based on a parameters
	a_order <- order(bank_info$PARAM.1)
	idx_astr <- as.numeric(cut(a_order, breaks=n.astr, labels=1:n.astr))
	ni_Va <- as.numeric(table(idx_astr))

	# assign item to each category of a-stratification
	Va <- list()
	for(i in 1:n.astr) {
		Va[[i]] <- c(1:I)[idx_astr == i]
	}

theta.mat <- array(NA, c(50, n.test))
se.mat <- array(NA, c(50, n.test))
item.mat <- array(NA, c(50, n.test))

for(r in 1:50) {

# make a model
shadow <- make.lp(nrow=0, ncol=I)

# set control parameters: maximization problem; integer tolerance is et to 0.1; 
# absolute MIP gaps is set to 0.1, relative MIP gap is set to 0.5
lp.control(lprec=shadow, sense="max", timeout=lp.control$timeout, epsint=lp.control$epsint, mip.gap=lp.control$mip.gap)

# contratint: type of decision variables
set.type(lprec=shadow, columns=1:I, type="binary")
set.bounds(lprec=shadow, lower=rep(0, I), upper=rep(1, I))

# constraint: test length
add.constraint(lprec=shadow, xt=rep(1, I), type="=", rhs=n.test, indices=1:I)

# constraints: content
for(i in 1:length(Nc)) {
	add.constraint(lprec=shadow, xt=rep(1, length(Vc[[i]])), type=">=", rhs=Nc[i], indices=Vc[[i]])
}

# constraints: a-stratification
for(i in 1:length(Na)) {
	add.constraint(lprec=shadow, xt=rep(1, length(Va[[i]])), type=">=", rhs=Na[i], indices=Va[[i]])
}


# 
theta.est <- ini.theta
item.recod <- c(); item.ans <- c(); theta.recod <- c(); se.recod <- c()

for(i in 1:n.test) {

	info <- as.numeric(test.info(x=item.pool, theta=theta.est, D=D)$itemInfo)
	
	# objective function
	set.objfn(lprec=shadow, obj=info)
	
	# constraint: include administered items
	if(i > 1) {
		add.constraint(lprec=shadow, xt=1, type="=", rhs=1, indices=item.recod[i - 1])
	}

	# solve teh model
	res_flag <- solve(shadow)
	
	# retrieve the values of the decision variables 
	sim_opt <- get.variables(shadow)

	# 
	item.admTotal <- c(1:I)[sim_opt == 1]
	Level <- as.numeric(cut(1:n.test, breaks=n.astr, level=1:n.astr))
	item.adm <- item.admTotal[idx_astr[sim_opt == 1] == Level[i]]
	
	if(i > 1) item.adm <- item.adm[!item.adm %in% item.recod]
	
	info.adm <- info[item.adm]
	if(NaCum[Level[i]] - (i - 1) >= ransq[Level[i]]) {
		item.ran <- item.adm[rank(info.adm) %in% 1:ransq[Level[i]]]
	} else {
		item.ran <- item.adm 
	}
	
	item.next <- ifelse(length(item.ran) > 1, sample(item.ran, 1), item.ran)
	item.recod <- c(item.recod, item.next)
	
	answer <- simdat(theta=true.theta, a.dc=a[item.next], b.dc=b[item.next], g.dc=g[item.next], cats=2, D=D)
	item.ans <- c(item.ans, answer)
	
	# intrim ability estimate
	my <- data.matrix(item.pool[item.recod, c("PARAM.1", "PARAM.2", "PARAM.3")])
	it <- cbind(my, rep(1, length(item.ans)))
	theta.est <- thetaEst(it=it, x=item.ans, D=D, method="EAP", priorDist="norm",
							priorPar=c(0, 1), range = c(-4, 4), parInt = c(-4, 4, 81))
	se.est <- semTheta(thEst=theta.est, it=it, x=item.ans, D=D, method = "EAP",
						priorDist = "norm", priorPar = c(0, 1), parInt = c(-4, 4, 81))
	theta.recod <- c(theta.recod, theta.est)
	se.recod <- c(se.recod, se.est)
	
}

theta.mat[r, ] <- theta.recod
se.mat[r, ] <- se.recod
item.mat[r, ] <- item.recod

print(r)

}




