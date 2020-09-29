
#####################################
##########  SET WORKSPACE  ##########
#####################################

rm(list = ls())  # fresh workspace

#####  Load Libraries  #####

# spatial packages
library(fields)  # rdist.earth function calculates great-circle distance
library(maps)  #

# parallelization packages
library(foreach)  # foreach function w/ %dopar% allows easy parallelization
library(doMC)  # used to register number of available cores

# graphics packages
library(ggplot2)
library(RColorBrewer)

#####  Load Functions  #####

source("functions.R")  # user-defined functions

#####  Load Data  #####

reformat <- function(A, make.matrix = FALSE) {
  rownames(A) <- A[, 1]  # set first column as rownames
  A <- A[, -1]  # remove first column
  if (make.matrix) A <- as.matrix(A)
  return(A)
}

S <- reformat(read.csv("data/S.csv", header = TRUE), make.matrix = TRUE)
X <- reformat(read.csv("data/X.csv", header = TRUE))
if (file.exists("data/Y.csv")) {
  Y <- reformat(read.csv("data/Y.csv", header = TRUE), make.matrix = TRUE)
} else {
  warning("File `data/Y.csv` not found, using `data/Y-reduced.csv` instead. This is a smaller version of the OTU table with only 5k taxa rather the original 50k+. The code will run faster, but the model results will not match the results from the full data.")
  Y <- reformat(read.csv("data/Y-reduced.csv", header = TRUE), make.matrix = TRUE)
}
if (file.exists("data/tax.csv")) {
  tax <- reformat(read.csv("data/tax.csv", header = TRUE), make.matrix = TRUE)
} else {
  warning("File `data/tax.csv` not found, using `data/tax-reduced.csv` instead. This is a smaller version of the taxonomy table, with only 5k taxa rather the original 50k+. The code will run faster, but the model results will not match the results from the full data.")
  tax <- reformat(read.csv("data/tax-reduced.csv", header = TRUE), make.matrix = TRUE)
}

rm(reformat)

#####  Session Info  #####

sessionInfo()

# R version 3.2.0 (2015-04-16)
# Platform: x86_64-apple-darwin13.4.0 (64-bit)
# Running under: OS X 10.11.2 (unknown)
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  grid      stats     graphics  grDevices utils     datasets  methods  
# [9] base     
# 
# other attached packages:
#   [1] RColorBrewer_1.1-2 ggplot2_2.0.0      doMC_1.3.3         iterators_1.0.7   
# [5] foreach_1.4.2      fields_8.3-5       maps_2.3-9         spam_1.0-1        
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_0.12.1      codetools_0.2-11 plyr_1.8.3.9000  gtable_0.1.2    
# [5] scales_0.3.0     tools_3.2.0      munsell_0.4.2    colorspace_1.2-6  