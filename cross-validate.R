
#####################################
#########  CROSS VALIDATE  ##########
#####################################

# make sure you have run source("get-data.R") at least once!
source("set-workspace.R")  # takes some time to load Y

########## Prepare to run model ##########

# Dimensions
n <- nrow(Y)
m <- ncol(Y)

# Specify number of cores to parallelize across
# It's best to run this file across many cores (>5)
ncore <- 10
registerDoMC(ncore)

# Additionally, split list of taxa by number of cores
taxa.list <- split(1:m, gl(ncore, ceiling(m / ncore), length = m))

# Cross-validation across nfold folds
nfold <- 5  # typically 5 or 10
fold <- rep(1:nfold, n)[1:n]
set.seed(777)
fold <- sample(fold)

# Possible kernel bandwidths (defined by rho in the accompanying manuscript)
bw <- c(100, 200, 300, 400, 500, 1000, 10000)

# Create prediction grid T across USA
# Finer resolution grid requires increased computation time
Tgrid <- generate_grid(x = 100, y = 100)

# Define empty list to store model results
results <- vector("list", nfold)

# record time before beginning cross-validation
ptm <- proc.time()

for (f in 1:nfold) {  # Begin cross-validation
  print(paste("Beginning fold", f, "of", nfold))
  
  # For fold f, split data into 80% train, 20% test
  # train used for two major steps (see below)
  # test will help determine performance of our model
  test   <- fold == f
  train  <- !test
  Ytest  <- Y[test, ]
  Stest  <- S[test, ]
  Ytrain <- Y[train, ]
  Strain <- S[train, ]
  
  # Split training data into 80% train2, 20% test2
  # train2 used for kernel smoothing & GCV
  # test2 used to find q for prediction regions
  set.seed(0820 * f)
  test2   <- runif(sum(train)) < 0.2
  train2  <- !test2
  Ytest2  <- Ytrain[test2, ]
  Stest2  <- Strain[test2, ]
  Ytrain2 <- Ytrain[train2, ]
  Strain2 <- Strain[train2, ]
  
  # Split Ytrain2, Ytest and Ytest2 by taxa.list for parallelization
  Ytrain2.list <- lapply(taxa.list, subset_columns, Ytrain2)
  Ytest.list <- lapply(taxa.list, subset_columns, Ytest)
  Ytest2.list <- lapply(taxa.list, subset_columns, Ytest2)
  
  # Train model via best bw and kernel smooth
  # Produces matrix of estimated occurence probabilities
  print("Estimating occurence probabilities...")
  M.list <- foreach(Y = Ytrain2.list) %dopar% 
              train_model(Y, S = Strain2, bw = bw, Tgrid = Tgrid)
  rm(Ytrain2.list, Strain2)
  
  ## Find log-likelihood for every sample in test and test2
  print("Calculating log-likelihood values for test set...")
  llike.test <- foreach(Y = Ytest.list, M = M.list, .combine = "+") %dopar% 
                  calculate_log_likelihood(Y, M)
  rm(Ytest.list)
  print("Calculating log-likelihood values for test2 set...")
  llike.test2 <- foreach(Y = Ytest2.list, M = M.list, .combine = "+") %dopar% 
                  calculate_log_likelihood(Y, M)
  rm(M.list, Ytest2.list)
     
  ## Obtain predictive pmfs by normalizing log-likelihoods
  pmf.test  <- make_pmf(llike.test)
  pmf.test2 <- make_pmf(llike.test2)
  rm(llike.test, llike.test2)
  
  ## Make spatial source predictions!
  Stest.hat <- Tgrid[apply(pmf.test, 2, which.max), ]
  rownames(Stest.hat) <- rownames(Stest)
  Stest2.hat <- Tgrid[apply(pmf.test2, 2, which.max), ]
  rownames(Stest2.hat) <- rownames(Stest2)

  ## Keep pmfs, test sets, and train sets
  print("Storing results...")
  results[[f]] <- list(pmf.test = pmf.test, Stest = Stest, Stest.hat = Stest.hat,
                       pmf.test2 = pmf.test2, Stest2 = Stest2, Stest2.hat = Stest2.hat)
  rm(pmf.test, Stest, Stest.hat, pmf.test2, Stest2, Stest2.hat)
  print(paste("Finished with fold", f, "of", nfold))
}  

# End cross-validation
proc.time() - ptm  # elapsed time of cross-validation

# save results
print("Saving final results...")
save(list = c("results", "Tgrid"), file = "results.RData")
print("Done!")
