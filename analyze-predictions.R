
#####################################
#######  ANALYZE PREDICTIONS  #######
#####################################

# make sure you have run source("get-data.R") at least once!
source("set-workspace.R")  # takes some time to load Y

if (!file.exists("results.RData")) {
  # Download results.RData locally if not in directory.
  # File size is 92.8 MB
  download.file("http://www4.ncsu.edu/~ngranth/results.RData", "results.RData", mode="wb")
}
load("results.RData")  # results from cross-validate.R

# So... How did we do?

# Recall that the errors here will not match the manuscript because 
# geographic coordinates have been truncated to protect the confidentiality
# of study participants. Thus, error is expected to be slightly higher than
# originally reported, even despite the inclusion of additional data
# (there are now 1,331 homes in our database vs the n = 928 in the manuscript)

pmf.test <- do.call(cbind, lapply(results, extract, "pmf.test"))
Stest <- do.call(rbind, lapply(results, extract, "Stest"))
Stest.hat <- do.call(rbind, lapply(results, extract, "Stest.hat"))

reorder.homes <- order(as.numeric(rownames(Stest)))
pmf.test <- pmf.test[, reorder.homes]
Stest <- Stest[reorder.homes, ]
Stest.hat <- Stest.hat[reorder.homes, ]
rm(reorder.homes)

## Let's investigate overall prediction error
pred.error <- diag(rdist.earth(Stest, Stest.hat, miles = FALSE))
summary(pred.error)
library(ggplot2)
p <- ggplot(as.data.frame(pred.error), aes(x = pred.error))
p <- p + geom_histogram(binwidth = 100)
p <- p + scale_x_continuous(breaks = seq(0, 4500, by = 500))
p <- p + xlab("Prediction Error (km)") + ylab("Count")
p  # Figure 4 in the manuscript

# Another way to judge the accuracy of our predictions is by how well 
# our prediction regions cover the true origin.

# Let a prediction region be defined by R_q <- {t in Tgrid: f(t) > k}, 
# where f(t) denotes the value of the predictive pmf at t and 
# k is a constant such that the sum of f(t) over t in R_q yields (at least) q.  

######  Find q using the test2 set (a subset of the train set)  #####

pmf.test2 <- lapply(results, extract, "pmf.test2")
Stest2 <- lapply(results, extract, "Stest2")
Stest2.hat <- lapply(results, extract, "Stest2.hat")

## Use test2 in each fold to find q for varying region probabilities
regions <- c(0.5, 0.75, 0.9)
# Takes several minutes to find optimal q in each fold for fine 'by'...
q.by.fold <- do.call(rbind, Map(select_q, pmf.test2, Stest2, 
                                MoreArgs = list(Tgrid = Tgrid, 
                                                regions = regions, 
                                                by = 0.005)))
q <- colMeans(q.by.fold)  # average over per-fold q for each region

######  Now use q to form pred regions in test set  #####

# How well do the q obtained on our test2 set retain coverage on test set?
calculate_coverage(pmf.test, q, Stest, Tgrid)  # fairly well

# What do our predictions look like?
rownames(Stest)  # valid homeIDs
plot_prediction("1227", pmf.test, q, Stest, Stest.hat, Tgrid, save = TRUE)  # Figure 3
# (worse performance here likely due to truncated geo coordinates and fold allocation of home 1227)
plot_prediction("44", pmf.test, q, Stest, Stest.hat, Tgrid, save = TRUE)
plot_prediction(c("515", "881"), pmf.test, q, Stest, Stest.hat, Tgrid, save = TRUE)
plot_prediction(c("480", "331", "173"), pmf.test, q, Stest, Stest.hat, Tgrid, save = TRUE)

## To generate ALL n prediction plots and save them in figs/predictions
## (size for all n prediction plots = ~600 Mb, runtime = ~54 minutes)
## uncomment the following code and run it:
# plot_prediction(rownames(Stest), pmf.test, q, Stest, Stest.hat, Tgrid, save = TRUE)

# It is also helpful to look at all predictions simultaneously to
# identify regional biases in our model.
plot_all_predictions(Stest, Stest.hat, save = TRUE)  # S1 Figure (Supplementary Information)

# Finally, how do the errors behave by covariates?
df.err <- cbind.data.frame(pred.error, is_covered(pmf.test, q, Stest, Tgrid))
summarize_error(df.err)  # overall

# By sampling intensity (i.e., number of neighbors within great-circle distance of 100 km)
num.neighbors <- rowSums(rdist.earth(Stest, miles = FALSE) <= 100) - 1
summary(num.neighbors)
hist(num.neighbors, breaks = 50)
summarize_error_by_variable(df.err, num.neighbors)

# By fungal richness (# of distinct taxa observed in a sample)
fungal.richness <- rowSums(Y)
summary(fungal.richness)
hist(fungal.richness, breaks = 50)
summarize_error_by_variable(df.err, fungal.richness)

# By mean temperature (above or below median)
mean.temp <- X$MeanAnnual.Temperature
summary(mean.temp)
hist(mean.temp, breaks = 50)
summarize_error_by_variable(df.err, mean.temp, probs = 0.5)

# By mean precipitation (above or below median)
mean.precip <- X$MeanAnnualPrecipitation
summary(mean.precip)
hist(mean.precip, breaks = 50)
summarize_error_by_variable(df.err, mean.precip, probs = 0.5)

# By elevation (above or below median)
elevation <- X$Elevation
summary(elevation)
hist(elevation, breaks = 50)
summarize_error_by_variable(df.err, elevation, probs = 0.5)
