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

##----------------------------------------------------------------------------
# Read all evaluation results (conditional results) 
eff_mfi <- readRDS(file.path(dir4save, "eff_mfi.rds"))
acc_mfi <- readRDS(file.path(dir4save, "acc_mfi.rds"))
eff_astr <- readRDS(file.path(dir4save, "eff_astr.rds"))
acc_astr <- readRDS(file.path(dir4save, "acc_astr.rds"))
eff_m1 <- readRDS(file.path(dir4save, "eff_m1.rds"))
acc_m1 <- readRDS(file.path(dir4save, "acc_m1.rds"))
eff_m2 <- readRDS(file.path(dir4save, "eff_m2.rds"))
acc_m2 <- readRDS(file.path(dir4save, "acc_m2.rds"))
eff_m3 <- readRDS(file.path(dir4save, "eff_m3.rds"))
acc_m3 <- readRDS(file.path(dir4save, "acc_m3.rds"))

# creat a data.frame to draw a plot
eff_df <- rbind(eff_mfi, eff_astr, eff_m1, eff_m2, eff_m2)
acc_df <- rbind(acc_mfi, acc_astr, acc_m1, acc_m2, acc_m2)


ggplot(my.df, aes(x = time, y = value) ) + 
  geom_line( aes(color = variable) ) + 
  facet_wrap(~Unit, scales = "free_y", nrow = 2, 
             strip.position = "left", 
             labeller = as_labeller(c(A = "Currents (A)", V = "Voltage (V)") ) )  +
  ylab(NULL) +
  theme(strip.background = element_blank(),
        strip.placement = "outside")



# This function draws plots
plot_shadow <- function(x, x.var="theta", y.var=c("CSEE", "Bias"), point.size=1.5, point.shape, line.size=0.8, line.color,
                        lab.size=15, axis.size=15, legend.size=15, strip.size=10, ylim=NULL, legend.title, 
                        legend.text, legend.position = "right", legend.key.width=1.5, legend.color.size=3) 

  

  


# This function draws plots for CSEE or Bias under Study 1
plot_study1 <- function(x, x.var="theta", y.var=c("CSEE", "Bias"), point.size=1.5, point.shape, line.size=0.8, line.color,
                        lab.size=15, axis.size=15, legend.size=15, strip.size=10, ylim=NULL, legend.title, 
                        legend.text, legend.position = "right", legend.key.width=1.5, legend.color.size=3) {
  
  if(y.var == "CSEE") ylab <- "Standard Error"
  if(y.var == "Bias") ylab <- "Bias"
  if(missing(line.color)) {
    line.color <-  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    # line.color <- brewer.pal(n = 9, name = "Greys")[3:9]
  } else {
    line.color <- line.color
  }
  if(missing(point.shape)) {
    point.shape <-  c(1, 2, 4, 3)
  } else {
    point.shape <- point.shape
  }
  if(missing(legend.title)) legend.title <- "Method"
  if(missing(legend.text)) legend.text <- levels(x$Method)
  
  if(is.null(ylim)) {
    p <- x %>% 
      ggplot(mapping=aes_string(x="theta", y=y.var)) +
      geom_point(mapping=aes_string(shape="Method"), size=point.size) +
      geom_line(mapping=aes_string(color="Method"), size=line.size) +
      # geom_line(mapping=aes_string(linetype="Method"), size=0.8) + 
      labs(x = expression(theta), y = ylab) +
      theme_bw() +
      facet_grid(Length ~ Panel) +
      theme(axis.title = element_text(size=lab.size),
            axis.text = element_text(size=axis.size)) +
      theme(legend.title = element_text(size=legend.size),
            legend.text = element_text(size=legend.size), 
            legend.position = legend.position) +
      theme(strip.text.x = element_text(size = strip.size, face = 'bold'),
            strip.text.y = element_text(size = strip.size, face = 'bold')) +
      scale_colour_manual(values=line.color, name = legend.title, labels = legend.text) +
      scale_shape_manual(values=point.shape, name = legend.title, labels = legend.text) +
      theme(legend.key.width = unit(1.5, "cm")) +
      guides(color = guide_legend(override.aes = list(size=legend.color.size)))
  } else {
    p <- x %>% 
      ggplot(mapping=aes_string(x="theta", y=y.var)) +
      geom_line(mapping=aes_string(color="Method"), size=line.size) +
      geom_point(mapping=aes_string(shape="Method"), size=point.size) +
      # geom_line(mapping=aes_string(linetype="Method"), size=0.8) + 
      labs(x = expression(theta), y = ylab) +
      ylim(ylim[1], ylim[2]) +
      theme_bw() +
      facet_grid(Length ~ Panel) +
      theme(axis.title = element_text(size=lab.size),
            axis.text = element_text(size=axis.size)) +
      theme(legend.title = element_text(size=legend.size),
            legend.text = element_text(size=legend.size), 
            legend.position = legend.position) +
      theme(strip.text.x = element_text(size = strip.size, face = 'bold'),
            strip.text.y = element_text(size = strip.size, face = 'bold')) +
      scale_colour_manual(values=line.color, name = legend.title, labels = legend.text) +
      scale_shape_manual(values=point.shape, name = legend.title, labels = legend.text) +
      theme(legend.key.width = unit(legend.key.width, "cm")) +
      guides(color = guide_legend(override.aes = list(size=legend.color.size)))
  }
  
  p
  
}

