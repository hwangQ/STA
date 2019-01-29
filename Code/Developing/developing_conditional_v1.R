###############################################################################
#Authors: Hwanggyu Lim (UMASS AMHERST)
#Date: 11/26/2018
#Title: Bottom-up assembly code 
###############################################################################
# call packages
library(lpSolveAPI)
library(reshape2)
library(parallel)
library(pbapply)
library(tidyverse)

##----------------------------------------------------------------------------
# Read source code
src.files <- list.files("R/", pattern="*.R")
purrr::walk(src.files, function(x) source(file.path("R", x)))

# read item bank and content category information
bank_info <- readRDS("Temp/item_bank_sim_nocor.rds")
content <- read.csv("Temp/Content_Requirement.csv")
colnames(content) <- c("CLASS", "MIN", "MAX")

# set condtions
nrep <- 500
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
theta.points <- seq(-3, 3, 0.5)
true.thetas <- purrr::map(theta.points, .f=function(x) rep(x, nrep))
ini.thetas <- 0
parallel <- TRUE
useSE <- FALSE
sd.con <- 0.5
source.files <- file.path(getwd(), "R", src.files)

################################################################################################
dir4save <- file.path("SaveFiles/Conditional")

##----------------------------------------------------------------------------------------------------------------
## 1. Maximum Fisher's Information method (baseline approach)
# (1) simulation
shadow_mfi <- ShadowCATall_parel(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas, ini.thetas=ini.thetas, 
                                 min.test=min.test, max.test=max.test, b.content=b.content, parallel=FALSE, useSE=FALSE, sd.con=sd.con,
                                 a.strat=FALSE, n.strat=n.strat, breaks=breaks, ransq=FALSE, par.ransq=par.ransq, post.ransq=post.ransq, 
                                 estimation=estimation, D=D, lp.control=lp.control, source.files, seed=107)

saveRDS(shadow_mfi, file.path(dir4save, "obj_shadow_mfi.rds"))
shadow_mfi <- readRDS(file.path(dir4save, "obj_shadow_mfi.rds"))

# (2) evaluation
# a function to evalute item utilization and item exposure rate
expose_mfi <- purrr::map(shadow_mfi, .f=exposRate, item.pool, crit.over=0.3, crit.less=0.02)

# a function to evalute efficiency of CAT
eff_mfi <- t(purrr::map_dfc(shadow_mfi, .f=function(x) cat_eff(x)[1:2])) %>% 
  data.frame() %>% 
  rename("mean.length"=X1, "unreached.rate"=X2) %>% 
  mutate(theta = theta.points, method = "MFI") %>% 
  reshape2::melt(variable.name="type", id.vars=c("theta", "method"), value.name="value")

# a function to evalute mesurement accuracy
acc_mfi <- t(purrr::map_dfc(shadow_mfi, .f=cat_acc)) %>% 
  data.frame() %>% 
  rename("bias"=X1, "rmse"=X2) %>% 
  mutate(theta = theta.points, method = "MFI") %>% 
  reshape2::melt(variable.name="type", id.vars=c("theta", "method"), value.name="value")

# a function to evalute the average of se
meanSE_mfi <- purrr::map_dfc(shadow_mfi, .f=function(x) mean(x$summary$final_df$SeEst))

# save evalution results
saveRDS(eff_mfi, file.path(dir4save, "eff_mfi.rds"))
saveRDS(acc_mfi, file.path(dir4save, "acc_mfi.rds"))
# eff_mfi <- readRDS(file.path(dir4save, "eff_mfi.rds"))
# acc_mfi <- readRDS(file.path(dir4save, "acc_mfi.rds"))

rm(shadow_mfi)


##----------------------------------------------------------------------------------------------------------------
## 2. a-stratification method (Diao & Ren (2018) approach)
# (1) simulation
shadow_astr <- ShadowCATall_parel(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas, ini.thetas=ini.thetas, 
                                 min.test=min.test, max.test=max.test, b.content=b.content, parallel=FALSE, useSE=FALSE, sd.con=sd.con,
                                 a.strat=TRUE, n.strat=n.strat, breaks=breaks, ransq=FALSE, par.ransq=par.ransq, post.ransq=post.ransq, 
                                 estimation=estimation, D=D, lp.control=lp.control, source.files, seed=172)

saveRDS(shadow_astr, file.path(dir4save, "obj_shadow_astr.rds"))
shadow_astr <- readRDS(file.path(dir4save, "obj_shadow_astr.rds"))

# (2) evaluation
# a function to evalute item utilization and item exposure rate
expose_astr <- purrr::map(shadow_astr, .f=exposRate, item.pool, crit.over=0.3, crit.less=0.02)

# a function to evalute efficiency of CAT
eff_astr <- t(purrr::map_dfc(shadow_astr, .f=function(x) cat_eff(x)[1:2])) %>% 
  data.frame() %>% 
  rename("mean.length"=X1, "unreached.rate"=X2) %>% 
  mutate(theta = theta.points, method = "ASTR") %>% 
  reshape2::melt(variable.name="type", id.vars=c("theta", "method"), value.name="value")

# a function to evalute mesurement accuracy
acc_astr <- t(purrr::map_dfc(shadow_astr, .f=cat_acc)) %>% 
  data.frame() %>% 
  rename("bias"=X1, "rmse"=X2) %>% 
  mutate(theta = theta.points, method = "ASTR") %>% 
  reshape2::melt(variable.name="type", id.vars=c("theta", "method"), value.name="value")

# a function to evalute the average of se
meanSE_astr <- purrr::map_dfc(shadow_astr, .f=function(x) mean(x$summary$final_df$SeEst))

# save evalution results
saveRDS(eff_astr, file.path(dir4save, "eff_astr.rds"))
saveRDS(acc_astr, file.path(dir4save, "acc_astr.rds"))
# eff_astr <- readRDS(file.path(dir4save, "eff_astr.rds"))
# acc_astr <- readRDS(file.path(dir4save, "acc_astr.rds"))

rm(shadow_astr)


##----------------------------------------------------------------------------------------------------------------
## 3. modification 1 
# (1) simulation
shadow_m1 <- ShadowCATall_parel(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas, ini.thetas=ini.thetas, 
                                  min.test=min.test, max.test=max.test, b.content=b.content, parallel=FALSE, useSE=FALSE, sd.con=sd.con,
                                  a.strat=TRUE, n.strat=n.strat, breaks=breaks, ransq=TRUE, par.ransq=par.ransq, post.ransq=post.ransq, 
                                  estimation=estimation, D=D, lp.control=lp.control, source.files, seed=329)

saveRDS(shadow_m1, file.path(dir4save, "obj_shadow_m1.rds"))
shadow_m1 <- readRDS(file.path(dir4save, "obj_shadow_m1.rds"))

# (2) evaluation
# a function to evalute item utilization and item exposure rate
expose_m1 <- purrr::map(shadow_m1, .f=exposRate, item.pool, crit.over=0.3, crit.less=0.02)

# a function to evalute efficiency of CAT
eff_m1 <- t(purrr::map_dfc(shadow_m1, .f=function(x) cat_eff(x)[1:2])) %>% 
  data.frame() %>% 
  rename("mean.length"=X1, "unreached.rate"=X2) %>% 
  mutate(theta = theta.points, method = "MODIF1") %>% 
  reshape2::melt(variable.name="type", id.vars=c("theta", "method"), value.name="value")

# a function to evalute mesurement accuracy
acc_m1 <- t(purrr::map_dfc(shadow_m1, .f=cat_acc)) %>% 
  data.frame() %>% 
  rename("bias"=X1, "rmse"=X2) %>% 
  mutate(theta = theta.points, method = "MODIF1") %>% 
  reshape2::melt(variable.name="type", id.vars=c("theta", "method"), value.name="value")

# a function to evalute the average of se
meanSE_m1 <- purrr::map_dfc(shadow_m1, .f=function(x) mean(x$summary$final_df$SeEst))

# save evalution results
saveRDS(eff_m1, file.path(dir4save, "eff_m1.rds"))
saveRDS(acc_m1, file.path(dir4save, "acc_m1.rds"))
# eff_m1 <- readRDS(file.path(dir4save, "eff_m1.rds"))
# acc_m1 <- readRDS(file.path(dir4save, "acc_m1.rds"))

rm(shadow_m1)


##----------------------------------------------------------------------------------------------------------------
## 4. modification 2 
# (1) simulation
shadow_m2 <- ShadowCATall_parel(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas, ini.thetas=ini.thetas, 
                                min.test=min.test, max.test=max.test, b.content=b.content, parallel=TRUE, useSE=FALSE, sd.con=sd.con,
                                a.strat=TRUE, n.strat=n.strat, breaks=breaks, ransq=TRUE, par.ransq=par.ransq, post.ransq=post.ransq, 
                                estimation=estimation, D=D, lp.control=lp.control, source.files, seed=772)

saveRDS(shadow_m2, file.path(dir4save, "obj_shadow_m2.rds"))
shadow_m2 <- readRDS(file.path(dir4save, "obj_shadow_m2.rds"))

# (2) evaluation
# a function to evalute item utilization and item exposure rate
expose_m2 <- purrr::map(shadow_m2, .f=exposRate, item.pool, crit.over=0.3, crit.less=0.02)

# a function to evalute efficiency of CAT
eff_m2 <- t(purrr::map_dfc(shadow_m2, .f=function(x) cat_eff(x)[1:2])) %>% 
  data.frame() %>% 
  rename("mean.length"=X1, "unreached.rate"=X2) %>% 
  mutate(theta = theta.points, method = "MODIF2") %>% 
  reshape2::melt(variable.name="type", id.vars=c("theta", "method"), value.name="value")

# a function to evalute mesurement accuracy
acc_m2 <- t(purrr::map_dfc(shadow_m2, .f=cat_acc)) %>% 
  data.frame() %>% 
  rename("bias"=X1, "rmse"=X2) %>% 
  mutate(theta = theta.points, method = "MODIF2") %>% 
  reshape2::melt(variable.name="type", id.vars=c("theta", "method"), value.name="value")

# a function to evalute the average of se
meanSE_m2 <- purrr::map_dfc(shadow_m2, .f=function(x) mean(x$summary$final_df$SeEst))

# save evalution results
saveRDS(eff_m2, file.path(dir4save, "eff_m2.rds"))
saveRDS(acc_m2, file.path(dir4save, "acc_m2.rds"))
# eff_m2 <- readRDS(file.path(dir4save, "eff_m2.rds"))
# acc_m2 <- readRDS(file.path(dir4save, "acc_m2.rds"))

rm(shadow_m2)


##----------------------------------------------------------------------------------------------------------------
## 5. modification 3 
# (1) simulation
shadow_m3 <- ShadowCATall_parel(item.pool=item.pool, crit.se=crit.se, true.thetas=true.thetas, ini.thetas=ini.thetas, 
                                min.test=min.test, max.test=max.test, b.content=b.content, parallel=TRUE, useSE=TRUE, sd.con=sd.con,
                                a.strat=TRUE, n.strat=n.strat, breaks=breaks, ransq=TRUE, par.ransq=par.ransq, post.ransq=post.ransq, 
                                estimation=estimation, D=D, lp.control=lp.control, seed=575)

saveRDS(shadow_m3, file.path(dir4save, "obj_shadow_m3.rds"))
shadow_m3 <- readRDS(file.path(dir4save, "obj_shadow_m3.rds"))

# (2) evaluation
# a function to evalute item utilization and item exposure rate
expose_m3 <- purrr::map(shadow_m3, .f=exposRate, item.pool, crit.over=0.3, crit.less=0.02)

# a function to evalute efficiency of CAT
eff_m3 <- t(purrr::map_dfc(shadow_m3, .f=function(x) cat_eff(x)[1:2])) %>% 
  data.frame() %>% 
  rename("mean.length"=X1, "unreached.rate"=X2) %>% 
  mutate(theta = theta.points, method = "MODIF3") %>% 
  reshape2::melt(variable.name="type", id.vars=c("theta", "method"), value.name="value")

# a function to evalute mesurement accuracy
acc_m3 <- t(purrr::map_dfc(shadow_m3, .f=cat_acc)) %>% 
  data.frame() %>% 
  rename("bias"=X1, "rmse"=X2) %>% 
  mutate(theta = theta.points, method = "MODIF3") %>% 
  reshape2::melt(variable.name="type", id.vars=c("theta", "method"), value.name="value")

# a function to evalute the average of se
meanSE_m3 <- purrr::map_dfc(shadow_m3, .f=function(x) mean(x$summary$final_df$SeEst))

# save evalution results
saveRDS(eff_m3, file.path(dir4save, "eff_m3.rds"))
saveRDS(acc_m3, file.path(dir4save, "acc_m3.rds"))
# eff_m3 <- readRDS(file.path(dir4save, "eff_m3.rds"))
# acc_m3 <- readRDS(file.path(dir4save, "acc_m3.rds"))

rm(shadow_m3)



