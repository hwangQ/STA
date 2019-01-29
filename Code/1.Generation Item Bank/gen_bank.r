# Read source code
src.dir <- "E:/Simulation/6.ShadowCAT/Rcode/R/"
src.files <- list.files(src.dir, pattern="*.r")
for(files in src.files) source(paste0(src.dir, files))

# set directory
mydir <- "E:/Simulation/6.ShadowCAT/Temp"
setwd(mydir)

## Set conditions for creating anchor item parameters
set.seed(123)
nitem <- 360
D <- 1

###################################################################################################
## 1. Bank: No correlation between a and b
# Create dichotomous item parameters
repeat{
	x <- rnorm(nitem, mean=1.0, sd=0.4)
	if(min(x) > 0.25 & max(x) < 2.2) {
		a <- x
		break
	} else {next}
}

repeat{
	x <- rnorm(nitem, mean=0, sd=1)
	if(min(x) > -3.0 & max(x) < 3.0) {
		b <- x
		break
	} else {next}
}

repeat{
	x <- rbeta(nitem, shape1=8, shape2=32)
	if(max(x) < 0.5) {
		g <- x
		break
	} else {next}
}

# content constraints using 4 categoies
categories <- c(rep("A", nitem/4), rep("B", nitem/4), rep("C", nitem/4), rep("D", nitem/4))
content <- sample(categories, size=nitem, replace = FALSE)
content_int <- as.numeric(as.factor(content))

# make a data frame of item bank
bank_info <- shape_df(par.dc=list(a=a, b=b, g=g), item.id=NULL, cats=2, model="3PLM")
bank_info$CLASS <- content
bank_info$CLASS_INT <- content_int

# save all parameter infomation
saveRDS(bank_info, "item_bank_sim_nocor.rds")
write.csv(bank_info, "item_bank_sim_nocor.csv")

###################################################################################################