
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
# (will be loaded in the pertinent plotting files
#  so as not to unnecessarily clutter the workspace)

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
Y <- reformat(read.csv("data/Y.csv", header = TRUE), make.matrix = TRUE)
tax <- reformat(read.csv("data/tax.csv", header = TRUE))

rm(reformat)

#####  Session Info  #####

sessionInfo()

# R version 3.1.2 (2014-10-31)
# Platform: x86_64-apple-darwin13.4.0 (64-bit)
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  grid      stats     graphics  grDevices utils     datasets  methods  
# [9] base     
# 
# other attached packages:
#   [1] doMC_1.3.3      iterators_1.0.7 fields_7.1      maps_2.3-7      spam_0.41-0    
# [6] foreach_1.4.2  
# 
# loaded via a namespace (and not attached):
#   [1] codetools_0.2-9 compiler_3.1.2  tools_3.1.2    