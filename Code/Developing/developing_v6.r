library(lpSolveAPI)
library(SuppDists)
library(catR)

##########################################################################
# Read source code
src.dir <- "E:/Simulation/6.ShadowCAT/Rcode/R/"
src.files <- list.files(src.dir, pattern="*.r")
for(files in src.files) source(paste0(src.dir, files))

# set directory
mydir <- "E:/Simulation/6.ShadowCAT/Temp"
setwd(mydir)

# read item bank and content category information
bank_info <- readRDS("item_bank_sim_nocor.rds")
content <- read.csv("Content_Requirement.csv")
colnames(content) <- c("CLASS", "MIN", "MAX")

# set condtions
item.pool <- bank_info
lp.control <- list(timeout=5, epsint=0.1, mip.gap=c(0.1, 0.05))
ini.theta <- 0
min.test <- 32
max.test <- 48
b.content <- list(lowb.content=content$MIN, upb.content=content$MAX)
D <- 1.702
a.strat <- TRUE
n.strat <- 3
breaks <- c(14, 10, 8)
ransq <- TRUE
par.ransq <- c(3, 2, 1)
post.ransq <- 2
inter <- list(method="WL", priorDist="norm", priorPar=c(0, 1), parInt = c(-4, 4, 81), range=c(-4, 4))
final <- list(method="ML", priorDist="norm", priorPar=c(0, 1), parInt = c(-4, 4, 81), range=c(-4, 4))
estimation <- list(inter=inter, final=final)
crit.se <- 0.2
set.seed(101)
true.thetas <- rnorm(50, 0, 1)
true.thetas <- ifelse(true.thetas > 3.5, 3.5, true.thetas)
true.thetas <- ifelse(true.thetas < -3.5, -3.5, true.thetas)
ini.thetas <- 0
parallel <- TRUE
useSE <- FALSE
sd.con <- 0.5
seed <- 321

##
x.0 <- ShadowCATall(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas, ini.thetas=ini.thetas, 
			min.test=min.test, max.test=max.test, b.content=b.content, parallel=parallel, useSE=useSE, sd.con=sd.con,
			a.strat=a.strat, n.strat=n.strat, breaks=breaks, ransq=ransq, par.ransq=par.ransq, post.ransq=post.ransq, 
			estimation=estimation, D=D, lp.control=lp.control, seed=seed)

			
			
x.1 <- ShadowCATall(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas, ini.thetas=ini.thetas, 
			min.test=min.test, max.test=max.test, b.content=b.content, parallel=parallel, useSE=useSE, sd.con=sd.con,
			a.strat=a.strat, n.strat=n.strat, breaks=breaks, ransq=ransq, par.ransq=par.ransq, post.ransq=1, 
			estimation=estimation, D=D, lp.control=lp.control, seed=seed)
			
x2 <- ShadowCATall(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas, ini.thetas=ini.thetas, 
			min.test=min.test, max.test=max.test, b.content=b.content, parallel=parallel, useSE=TRUE, sd.con=sd.con,
			a.strat=a.strat, n.strat=n.strat, breaks=breaks, ransq=ransq, par.ransq=par.ransq, post.ransq=post.ransq, 
			estimation=estimation, D=D, lp.control=lp.control, seed=seed)
			
			
y <- ShadowCATall(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas, ini.thetas=ini.thetas, 
			min.test=min.test, max.test=max.test, b.content=b.content, parallel=FALSE, useSE=FALSE, sd.con=sd.con,
			a.strat=TRUE, n.strat=n.strat, breaks=breaks, ransq=FALSE, par.ransq=par.ransq, post.ransq=1, 
			estimation=estimation, D=D, lp.control=lp.control, seed=seed)

## Evaluation
x <- x.0

# a function to evalute item utilization and item exposure rate
exposRate(x, item.pool, crit.over=0.3, crit.less=0.02)

# a function to evalute efficiency of CAT
cat_eff(x)

# a function to evalute mesurement accuracy
cat_acc(x)



