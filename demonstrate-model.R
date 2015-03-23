
#####################################
########  DEMONSTRATE MODEL  ########
#####################################

# make sure you have run source("get-data.R") at least once!
source("set-workspace.R")  # takes some time to load Y

# The purpose of this file is to showcase the model in action on a
# much smaller, more manageable data set. To this end, we restrict 
# the number of species to the first 1,000.

Y <- Y[, 1:1000]
n <- nrow(Y)
m <- ncol(Y)

########## Prepare to run model ##########

# Cross-validation across nfold folds
nfold <- 5  # typically 5 or 10
fold <- rep(1:nfold, n)[1:n]
set.seed(777)
fold <- sample(fold)
testing.group.by.fold <- lapply(1:nfold, function(f) fold == f)

# Kernel bandwidths (defined by rho in the accompanying manuscript)
bw <- c(100, 200, 300, 400, 500, 1000, 10000)

# Create prediction grid T across USA
# Finer resolution grid requires increased computation time
Tgrid <- generate_grid(x = 100, y = 100)

##########  Run the model!  ##########

f <- 1  # We will only run across the first fold in this file!

# For fold f, split data into 80% train, 20% test
# train used for two major steps (see below)
# test will help determine performance of our model
test   <- fold == f
train  <- !test
ntest  <- sum(test)
ntrain <- sum(train)
Ytest  <- Y[test, ]
Stest  <- S[test, ]
Ytrain <- Y[train, ]
Strain <- S[train, ]

# Split training data into 80% train2, 20% test2
# train2 used for kernel smoothing & GCV
# test2 used to find q for prediction regions
set.seed(0820)
test2   <- runif(sum(train)) < 0.2
train2  <- !test2
ntest2  <- sum(test2)
ntrain2 <- sum(train2)
Ytest2  <- Ytrain[test2, ]
Stest2  <- Strain[test2, ]
Ytrain2 <- Ytrain[train2, ]
Strain2 <- Strain[train2, ]

# Train model via best bw and kernel smooth
# Produces matrix of estimated occurence probabilities
M <- train_model(Y = Ytrain2, S = Strain2, bw = bw, Tgrid = Tgrid)

# Next, get log likelihood values for every sample in test and test2
llike.test <- calculate_log_likelihood(Ytest, M)
llike.test2 <- calculate_log_likelihood(Ytest2, M)

# Normalize log-likelihood values (necessary for forming prediction regions)
pmf.test <- make_pmf(llike.test)
pmf.test2 <- make_pmf(llike.test2)

# predict the origin that maximizes a sample's log-likelihood
Stest.hat <- Tgrid[apply(pmf.test, 2, which.max), ]
Stest2.hat <- Tgrid[apply(pmf.test2, 2, which.max), ]

# How did we do?
pred.error <- diag(rdist.earth(Stest, Stest.hat, miles = FALSE))
summary(pred.error)
hist(pred.error, breaks = 50)

# Not exceptional, but recall:
# 1. we restricted to only 1,000 species for the purposes of this example.
#    In reality we have >40k.
# 2. the locations in S are not exact and have been truncated to avoid
#    identifying the homes of participants in our study

# Further investigation of prediction error:
## Use cross-validation on test2 to find q for varying regions
regions <- c(0.5, 0.75, 0.9)
q <- select_q(pmf.test2, Stest2, Tgrid, regions, by = 0.005)

# How well do the q obtained on our test2 set retain coverage on test set?
calculate_coverage(pmf.test, q, Stest, Tgrid)  # we hope for 0.5, 0.75, 0.9

# What do our predictions look like?
plot_prediction(50, pmf.test, q, Stest, Stest.hat, Tgrid)
plot_prediction(100, pmf.test, q, Stest, Stest.hat, Tgrid)
plot_prediction(200, pmf.test, q, Stest, Stest.hat, Tgrid)
