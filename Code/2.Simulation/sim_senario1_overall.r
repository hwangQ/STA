################################################################################################
# Read source code
src.dir <- "E:/Simulation/6.ShadowCAT/Rcode/R/"
src.files <- list.files(src.dir, pattern="*.R")
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
breaks <- c(14, 12, 6)
ransq <- TRUE
par.ransq <- c(4, 3, 2)
post.ransq <- 1
inter <- list(method="WL", priorDist="norm", priorPar=c(0, 1), parInt = c(-4, 4, 81), range=c(-4, 4))
final <- list(method="ML", priorDist="norm", priorPar=c(0, 1), parInt = c(-4, 4, 81), range=c(-4, 4))
estimation <- list(inter=inter, final=final)
crit.se <- 0.2
set.seed(107)
true.thetas <- rnorm(3000, 0, 1)
true.thetas <- ifelse(true.thetas > 3.5, 3.5, true.thetas)
true.thetas <- ifelse(true.thetas < -3.5, -3.5, true.thetas)
ini.thetas <- 0
parallel <- TRUE
useSE <- FALSE
sd.con <- 0.5

################################################################################################
dir4save <- file.path("SaveFiles/Overall")

##----------------------------------------------------------------------------------------------------------------
## 1. Maximum Fisher's Information method (baseline approach)
# (1) simulation
shadow_mfi <- ShadowCATall(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas, ini.thetas=ini.thetas, 
					min.test=min.test, max.test=max.test, b.content=b.content, parallel=FALSE, useSE=FALSE, sd.con=sd.con,
					a.strat=FALSE, n.strat=n.strat, breaks=breaks, ransq=FALSE, par.ransq=par.ransq, post.ransq=post.ransq, 
					estimation=estimation, D=D, lp.control=lp.control, seed=273)

saveRDS(shadow_mfi, file.path(dir4save, "obj_shadow_mfi.rds"))
shadow_mfi <- readRDS(file.path(dir4save, "obj_shadow_mfi.rds"))

# (2) evaluation
# a function to evalute item utilization and item exposure rate
expose_mfi <- exposRate(shadow_mfi, item.pool, crit.over=0.3, crit.less=0.02)

# a function to evalute efficiency of CAT
eff_mfi <- cat_eff(shadow_mfi)

# a function to evalute mesurement accuracy
acc_mfi <- cat_acc(shadow_mfi)

# all of evaluation
eval_mfi <- unlist(c(acc_mfi, expose_mfi$over.expos[2], expose_mfi$less.expos[2], expose_mfi$unused[2], 
					expose_mfi$index.expos, expose_mfi$overlap.rate, eff_mfi))
names(eval_mfi) <- c("BIAS", "RMSE", "OVER.EXPOS", "LESS.EXPOS", "UNUSED", "CHISQ", 
					"OVERLAB", "MEAN.LENGTH", "UNREACHED", "EFFICIENCY") 

saveRDS(eval_mfi, file.path(dir4save, "eval_shadow_mfi.rds"))
rm(shadow_mfi)

##----------------------------------------------------------------------------------------------------------------
## 2. a-stratification method (Diao & Ren (2018) approach)
# (1) simulation
shadow_astr <- ShadowCATall(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas, ini.thetas=ini.thetas, 
					min.test=min.test, max.test=max.test, b.content=b.content, parallel=FALSE, useSE=FALSE, sd.con=sd.con,
					a.strat=TRUE, n.strat=n.strat, breaks=breaks, ransq=FALSE, par.ransq=par.ransq, post.ransq=post.ransq, 
					estimation=estimation, D=D, lp.control=lp.control, seed=357)

saveRDS(shadow_astr, file.path(dir4save, "obj_shadow_astr.rds"))
shadow_astr <- readRDS(file.path(dir4save, "obj_shadow_astr.rds"))

# (2) evaluation
# a function to evalute item utilization and item exposure rate
expose_astr <- exposRate(shadow_astr, item.pool, crit.over=0.3, crit.less=0.02)

# a function to evalute efficiency of CAT
eff_astr <- cat_eff(shadow_astr)

# a function to evalute mesurement accuracy
acc_astr <- cat_acc(shadow_astr)

# all of evaluation
eval_astr <- unlist(c(acc_astr, expose_astr$over.expos[2], expose_astr$less.expos[2], expose_astr$unused[2], 
					expose_astr$index.expos, expose_astr$overlap.rate, eff_astr))
names(eval_astr) <- c("BIAS", "RMSE", "OVER.EXPOS", "LESS.EXPOS", "UNUSED", "CHISQ", 
					"OVERLAB", "MEAN.LENGTH", "UNREACHED", "EFFICIENCY") 

saveRDS(eval_astr, file.path(dir4save, "eval_shadow_astr.rds"))
rm(shadow_astr)


##----------------------------------------------------------------------------------------------------------------
## 3. modification 1 
# (1) simulation
shadow_m1 <- ShadowCATall(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas, ini.thetas=ini.thetas, 
					min.test=min.test, max.test=max.test, b.content=b.content, parallel=FALSE, useSE=FALSE, sd.con=sd.con,
					a.strat=TRUE, n.strat=n.strat, breaks=breaks, ransq=TRUE, par.ransq=par.ransq, post.ransq=post.ransq, 
					estimation=estimation, D=D, lp.control=lp.control, seed=377)

saveRDS(shadow_m1, file.path(dir4save, "obj_shadow_m1.rds"))
shadow_m1 <- readRDS(file.path(dir4save, "obj_shadow_m1.rds"))

# (2) evaluation
# a function to evalute item utilization and item exposure rate
expose_m1 <- exposRate(shadow_m1, item.pool, crit.over=0.3, crit.less=0.02)

# a function to evalute efficiency of CAT
eff_m1 <- cat_eff(shadow_m1)

# a function to evalute mesurement accuracy
acc_m1 <- cat_acc(shadow_m1)

# all of evaluation
eval_m1 <- unlist(c(acc_m1, expose_m1$over.expos[2], expose_m1$less.expos[2], expose_m1$unused[2], 
					expose_m1$index.expos, expose_m1$overlap.rate, eff_m1))
names(eval_m1) <- c("BIAS", "RMSE", "OVER.EXPOS", "LESS.EXPOS", "UNUSED", "CHISQ", 
					"OVERLAB", "MEAN.LENGTH", "UNREACHED", "EFFICIENCY") 

saveRDS(eval_m1, file.path(dir4save, "eval_shadow_m1.rds"))
rm(shadow_m1)


##----------------------------------------------------------------------------------------------------------------
## 4. modification 2 
# (1) simulation
shadow_m2 <- ShadowCATall(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas, ini.thetas=ini.thetas, 
					min.test=min.test, max.test=max.test, b.content=b.content, parallel=TRUE, useSE=FALSE, sd.con=sd.con,
					a.strat=TRUE, n.strat=n.strat, breaks=breaks, ransq=TRUE, par.ransq=par.ransq, post.ransq=post.ransq, 
					estimation=estimation, D=D, lp.control=lp.control, seed=727)
					
saveRDS(shadow_m2, file.path(dir4save, "obj_shadow_m2.rds"))
shadow_m2 <- readRDS(file.path(dir4save, "obj_shadow_m2.rds"))

# (2) evaluation
# a function to evalute item utilization and item exposure rate
expose_m2 <- exposRate(shadow_m2, item.pool, crit.over=0.3, crit.less=0.02)

# a function to evalute efficiency of CAT
eff_m2 <- cat_eff(shadow_m2)

# a function to evalute mesurement accuracy
acc_m2 <- cat_acc(shadow_m2)

# all of evaluation
eval_m2 <- unlist(c(acc_m2, expose_m2$over.expos[2], expose_m2$less.expos[2], expose_m2$unused[2], 
					expose_m2$index.expos, expose_m2$overlap.rate, eff_m2))
names(eval_m2) <- c("BIAS", "RMSE", "OVER.EXPOS", "LESS.EXPOS", "UNUSED", "CHISQ", 
					"OVERLAB", "MEAN.LENGTH", "UNREACHED", "EFFICIENCY") 

saveRDS(eval_m2, file.path(dir4save, "eval_shadow_m2.rds"))
rm(shadow_m2)


##----------------------------------------------------------------------------------------------------------------
## 5. modification 3 
# (1) simulation
shadow_m3 <- ShadowCATall(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas, ini.thetas=ini.thetas, 
					min.test=min.test, max.test=max.test, b.content=b.content, parallel=TRUE, useSE=TRUE, sd.con=sd.con,
					a.strat=TRUE, n.strat=n.strat, breaks=breaks, ransq=TRUE, par.ransq=par.ransq, post.ransq=post.ransq, 
					estimation=estimation, D=D, lp.control=lp.control, seed=727)

saveRDS(shadow_m3, file.path(dir4save, "obj_shadow_m3.rds"))
shadow_m3 <- readRDS(file.path(dir4save, "obj_shadow_m3.rds"))

# (2) evaluation
# a function to evalute item utilization and item exposure rate
expose_m3 <- exposRate(shadow_m3, item.pool, crit.over=0.3, crit.less=0.02)

# a function to evalute efficiency of CAT
eff_m3 <- cat_eff(shadow_m3)

# a function to evalute mesurement accuracy
acc_m3 <- cat_acc(shadow_m3)

# all of evaluation
eval_m3 <- unlist(c(acc_m3, expose_m3$over.expos[2], expose_m3$less.expos[2], expose_m3$unused[2], 
					expose_m3$index.expos, expose_m3$overlap.rate, eff_m3))
names(eval_m3) <- c("BIAS", "RMSE", "OVER.EXPOS", "LESS.EXPOS", "UNUSED", "CHISQ", 
					"OVERLAB", "MEAN.LENGTH", "UNREACHED", "EFFICIENCY") 

saveRDS(eval_m3, file.path(dir4save, "eval_shadow_m3.rds"))
rm(shadow_m3)


##----------------------------------------------------------------------------------------------------------------
## calculate mean se
shadow_mfi <- readRDS(file.path(dir4save, "obj_shadow_mfi.rds"))
meanSE_mfi <- mean(shadow_mfi$summary$final_df$SeEst)
rm(shadow_mfi)

shadow_astr <- readRDS(file.path(dir4save, "obj_shadow_astr.rds"))
meanSE_astr <- mean(shadow_astr$summary$final_df$SeEst)
rm(shadow_astr)

shadow_m1 <- readRDS(file.path(dir4save, "obj_shadow_m1.rds"))
meanSE_m1 <- mean(shadow_m1$summary$final_df$SeEst)
rm(shadow_m1)

shadow_m2 <- readRDS(file.path(dir4save, "obj_shadow_m2.rds"))
meanSE_m2 <- mean(shadow_m2$summary$final_df$SeEst)
rm(shadow_m2)

shadow_m3 <- readRDS(file.path(dir4save, "obj_shadow_m3.rds"))
meanSE_m3 <- mean(shadow_m3$summary$final_df$SeEst)
rm(shadow_m3)

meanSE <- c(meanSE_mfi, meanSE_astr, meanSE_m1, meanSE_m2, meanSE_m3)

## summary table
sim1_table <- rbind(eval_mfi, eval_astr, eval_m1, eval_m2, eval_m3)
sim1_table <- round(cbind(sim1_table, meanSE), 3)

# write.csv(sim1_table, file.path(dir4save, "sim1_allEval_table.csv"))




