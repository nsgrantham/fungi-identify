
#####################################
###########  FUNCTIONS  #############
#####################################

########## General functions ##########

subset_columns <- function(cols, A) {
  stopifnot(all(cols %in% 1:ncol(A)))
  return(A[, cols])
}

find_row <- function(name, all.names) {
  row <- which(all.names %in% name)
  if (!any(row)) row <- NA
  return(row)
}

collapse_name <- function(x) {
  paste(x[!is.na(x)], collapse = "; ")
}

extract <- function(list, var) {
  list[[var]]
}

########## US Geography ##########

generate_grid <- function(x = 100, y = 100){
  # purpose: generate grid of points over the US at which to make predictions
  # returns: two column matrix, longitude and latitude of all grid points
  library(maps)  # for map.where
  corners <- enclose_USA()
  # create grid
  grid <- expand.grid(seq(corners[1], corners[2], length = x), 
                      seq(corners[3], corners[4], length = y))
  # retain only points that fall over US soil
  inUSA <- !is.na(map.where("usa", x = grid[, 1], y = grid[, 2]))
  grid <- as.matrix(grid[inUSA, ])
  # name the columns and rows
  colnames(grid) <- c("lon", "lat")
  rownames(grid) <- 1:nrow(grid)
  return(grid)
}

enclose_USA <- function() {
  # purpose: geographic coordinate limits within the continental US
  # returns: vector of corners of USA
  max.lat <- 49.384472    # most northern
  min.lat <- 24.520833    # most southern
  max.lon <- -66.947028   # most eastern
  min.lon <- -124.733056  # most western
  corners <- c(min.lon, max.lon, min.lat, max.lat)
}

########## Training & Running the Model ##########

get_weights <- function(S.row, S.col, bw) {  
  # purpose: calculate w_ij as given in eq. (1) in Supplementary Material for each bandwidth
  # returns: list of matrices of weights between S.row and S.col for each bandwidth bw
  library(fields)  # for rdist.earth
  dist <- rdist.earth(S.row, S.col, miles = FALSE)  # great-circle distance in km
  kern <- lapply(bw, function(b) exp(-0.5 * (dist / b)^2))  # gaussian
  w <- lapply(kern, function(K) K / rowSums(K))  # normalize
  return(w)
}

squeeze <- function(A, tol = 1e-4) {
  # purpose: bound values between tol and 1 - tol for computational stability
  # returns: A after low/high values are squeezed to tol/1-tol respectively
  A[A < tol]     <- tol
  A[A > 1 - tol] <- 1 - tol
  return(A)
}

select_best_bandwidth <- function(Y, S, bw) {
  # purpose: select the best bandwidth rho for use in the gaussian kernel smoother
  # returns: scalar, the best minimizer of GCV in bw
  Y <- as.matrix(Y)
  n <- nrow(Y)
  stopifnot(nrow(S) == n)
  w    <- get_weights(S, S, bw)  # matrix W_rho for each rho
  trc  <- lapply(w, function(w) sum(diag(w)))  # trace(W_rho) for each rho
  yhat <- lapply(w, function(w) squeeze(w %*% Y))  # yhat for each rho
  sse  <- lapply(yhat, function(y) colSums((Y - y)^2))
  gcv  <- do.call(cbind, Map(function(t, s) s / (n*(1 - t/n)^2), trc, sse))
  best <- apply(as.matrix(gcv), 1, which.min)  # choose bw that minimizes gcv(bw)
  return(best)
}

kernel_smoother <- function(Y, S, bw, best, Tgrid) {
  # purpose: apply gaussian kernel smoother to all values in Y over Tgrid
  # returns: matrix with estimated occurence probabilities for all N prediction
  #          locations (rows) and each species (column)
  Y <- as.matrix(Y)
  stopifnot(nrow(Y) == nrow(S))
  stopifnot(ncol(Y) == length(best))
  species.by.bw <- split(1:ncol(Y), factor(best, levels = 1:length(bw)))
  Y.by.bw <- lapply(species.by.bw, function(j) Y[, j])
  w <- get_weights(Tgrid, S, bw)
  M <- do.call(cbind, Map("%*%", w, Y.by.bw))[, order(unlist(species.by.bw))]
  return(squeeze(M))
}

train_model <- function(Y, S, bw, Tgrid) {
  # purpose: train the spatial prediction model
  # returns: matrix of estimated occurrence probabilties across Tgrid
  stopifnot(nrow(Y) == nrow(S))
  best <- select_best_bandwidth(Y, S, bw)
  M <- kernel_smoother(Y, S, bw, best, Tgrid)
  return(M)
}

calculate_log_likelihood <- function(Y, M) {
  # purpose: calculate log likelihood of many samples Y given est. occurrence probs M
  # returns: matrix of log likelihood values of Y 
  ll <- apply(Y, MARGIN = 1, calculate_log_likelihood_helper, M = M)
  return(ll)
}

calculate_log_likelihood_helper <- function(Y, M) {  
  # purpose: calculate log likelihood for single sample Y given est. occurrence probs M 
  # returns: log likelihood that locs in Tgrid produced Y
  Y.mat <- matrix(Y, nrow(M), length(Y), byrow = TRUE)
  ll    <- rowSums(dbinom(Y.mat, 1, M, log = TRUE))
  return(ll)
}

make_pmf <- function(ll) {
  # purpose: normalize matrix of log likelhood values
  # returns: normalized matrix of log likelihood values
  ll  <- as.matrix(ll)
  pos <- ll - min(ll)
  pmf <- t(t(pos) / colSums(pos))
  return(pmf)
}


########## Form Prediction Regions ##########

find_threshold <- function(pmf, q)  {
  # purpose: for pmf f, find k s.t. sum of f(t) over all t in Rq is at least q
  #          where Rq = {t in Tgrid : f(t) > k}
  # returns: k, the pmf threshold value
  pmf <- rev(sort(pmf/sum(pmf)))            # normalize, sort in dec. order
  k   <- pmf[min(which(cumsum(pmf) >= q))]  # sum of pmf geq c is at least q
  return(k)
}

is_in_Rq <- function(pmf, q) {
  k <- apply(pmf, 2, find_threshold, q)
  return(t(t(pmf) >= k))
}

is_closest_t <- function(S, Tgrid) {
  Tgrid.vs.S <- rdist.earth(Tgrid, S, miles = FALSE)
  return(apply(Tgrid.vs.S, 2, function(x) x == min(x)))
}

is_covered <- function(pmf, q, S, Tgrid) {
  t.bool  <- is_closest_t(S, Tgrid)  # locate t in Tgrid closest to each origin s
  # (treat this as the "true origin" for purposes of determining coverage)
  Rq.bool <- lapply(q, is_in_Rq, pmf = pmf)  # which values of the pmf exceed threshold for each q
  # finally, what is the proportion of closest t covered by Rq for each q?
  covered <- do.call(cbind, lapply(Rq.bool, function(bool) colSums(bool & t.bool)))
  return(covered)
}

calculate_coverage <- function(pmf, q, S, Tgrid) {
  covered <- is_covered(pmf, q, S, Tgrid)
  coverage <- colMeans(covered)
  return(coverage)
}

select_q <- function(pmf, S, Tgrid, regions, by) {
  stopifnot(nrow(pmf) == nrow(Tgrid))  # rows of pmf must be 1-to-1 with locs in Tgrid
  stopifnot(ncol(pmf) == nrow(S))  # cols of pmf must be 1-to-1 with locs in S
  q.possible <- seq(by, 1 - by, by = by)
  coverage <- calculate_coverage(pmf, q = q.possible, S = S, Tgrid = Tgrid)
  regions <- sort(regions)
  q.regions <- sapply(regions, function(prop) q.possible[min(which(coverage >= prop))])
  names(q.regions) <- regions
  return(q.regions)
}

form_prediction_regions <- function(pmf, q, Tgrid) {
  k <- sapply(q, find_threshold, pmf = pmf)
  in.regions <- t(t(matrix(rep(pmf, length(k)), ncol = length(k))) >= k)
  region <- rowSums(in.regions)
  Rq <- cbind.data.frame(Tgrid, factor(region))
  colnames(Rq) <- c("lon", "lat", "region")
  Rq <- Rq[Rq$region != 0, ]  # remove all locations not in any pred regions
  return(Rq)
}

########## Plotting Occurrence & Predictions ##########

# define function that maps occurence probabilities of a fungal OTU
plot_occurrence <- function(otu, tax, bw, Tgrid, save = FALSE, path = "figs/occurrence") {
  library(ggplot2)
  library(RColorBrewer)
  dir.create(path)
  .e <- environment()
  otu.match <- sapply(otu, find_row, rownames(tax))
  otu <- otu[!is.na(otu.match)]
  otu.match <- otu.match[!is.na(otu.match)]
  stopifnot(any(otu.match))  # stop if no valid matches
  otu.name <- apply(tax[otu.match, ], 1, collapse_name)
  print(paste("Generating", length(otu.match), "occurrence map(s)..."))
  M <- as.matrix(train_model(Y[, otu.match], S, bw, Tgrid))
  for (j in 1:ncol(M)) {
    shading <- cbind.data.frame(Tgrid, M[, j])
    colnames(shading) <- c("lon", "lat", "prob")
    samples <- cbind.data.frame(S, Y[, otu.match[j]])
    colnames(samples) <- c("lon", "lat", "present")
    color <- brewer.pal(9, "BuPu")[c(3, 9)]
    p <- ggplot(data = shading, aes(x = lon, y = lat, fill = prob), environment = .e)
    p <- p + geom_tile()
    p <- p + scale_fill_gradient(low = color[1], high = color[2], space = "Lab", 
                                 limits = c(0, 1), breaks = c(0.2, 0.4, 0.6, 0.8),
                                 name = "Occurence\nProbability")
    p <- p + geom_polygon(data = map_data("state"), aes(x = long, y = lat, group = group), 
                          colour = "black", fill = "white", alpha = 0)
    p <- p + geom_point(data = samples, aes(x = lon, y = lat, shape = factor(1 - present), alpha = factor(present)), 
                        size = 3, environment = .e, inherit.aes = FALSE)
    p <- p + scale_shape_manual(values = c(1, 4), name = "Present?", labels = c("Yes", "No"))
    p <- p + scale_alpha_manual(values = c(0.5, 1), guide = FALSE)
    p <- p + guides(shape = guide_legend(order = 1))
    p <- p + ggtitle(otu.name[j])
    p <- p + theme_bw()
    p <- p + theme(plot.title = element_text(size = 10), plot.background = element_blank(),
                   axis.line = element_blank(), axis.text.x = element_blank(),
                   axis.text.y = element_blank(), axis.ticks = element_blank(),
                   axis.title.x = element_blank(), axis.title.y = element_blank(),
                   panel.background = element_blank(), panel.border =element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(p) 
    if (save) ggsave(paste0(otu[j], ".png"), path = path)
  }
}

plot_prediction <- function(homeID, pmf, q, S, S.hat, Tgrid, save = FALSE, path = "figs/predictions") {
  library(ggplot2)
  library(RColorBrewer)
  dir.create(path)
  stopifnot(length(q) <= 9)  # number of distinct regions limited by color palette
  home.match <- sapply(homeID, find_row, rownames(S))
  homeID <- homeID[!is.na(home.match)]
  home.match <- home.match[!is.na(home.match)]
  stopifnot(any(home.match))  # stop if no valid matches
  print(paste("Generating", length(home.match), "prediction map(s)..."))
  rownames(S.hat) <- paste0(rownames(S), ".hat")
  for (i in 1:length(home.match)) {
    ## Predictions Regions
    Rq <- form_prediction_regions(pmf[, home.match[i]], q, Tgrid)
    p <- ggplot(data = Rq, aes(x = lon, y = lat, fill = region))
    p <- p + geom_tile()
    p <- p + scale_fill_manual(name = "Prediction\nRegion", 
                               labels = rev(names(q)), 
                               values = rev(rev(brewer.pal(9, "BuPu"))[1:length(q)]))
    ## Add US States
    p <- p + geom_polygon(data = map_data("state"), aes(x = long, y = lat, group = group), 
                          colour = "black", fill = "white", alpha = 0)
    ## True and Predicted Origin
    points <- as.data.frame(rbind(S[home.match[i], ], S.hat[home.match[i], ]))
    colnames(points) <- c("lon", "lat")
    p <- p + geom_point(data = points, aes(x = lon, y = lat, colour = factor(2:1)), 
                        size = 3, inherit.aes = FALSE)
    ## Connect with a line
    p <- p + geom_line(data = points, aes(x = lon, y = lat, group = factor(c(1, 1))), inherit.aes = FALSE)
    ## Modify legend & theme
    p <- p + scale_colour_discrete(name = "Origin", label = c("Predicted", "True"))  # labels
    p <- p + guides(colour = guide_legend(order = 1)) # put Origin above Probability
    p <- p + ggtitle(paste("Error:", round(drop(rdist.earth(points[1, ], points[2, ], miles = FALSE)), 1), "km"))
    p <- p + theme(axis.line = element_blank(), axis.text.x = element_blank(),
                   axis.text.y = element_blank(), axis.ticks = element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   panel.background = element_blank(),
                   panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.background = element_blank())
    print(p)  
    if (save) ggsave(paste0(homeID[i], ".png"), path = path)
  }
}

plot_all_predictions <- function(S, S.hat, save = FALSE, path = "figs") {
  stopifnot(nrow(S) == nrow(S.hat))
  library(ggplot2)
  n <- nrow(S)
  is.true.origin <- factor(c(rep(1, n), rep(0, n)))
  pairing <- factor(rep(1:n, 2))
  rownames(S.hat) <- paste0(rownames(S), ".hat")
  points <- cbind.data.frame(rbind(S, S.hat), is.true.origin, pairing)
  colnames(points) <- c("lon", "lat", "true", "pair")
  ## Plot points
  p <- ggplot(data = points, aes(x = lon, y = lat))
  p <- p + geom_point(aes(colour = true), size = 2, alpha = 0.8)
  ## Connect each true origin with its predicted origin
  p <- p + geom_line(aes(group = pair), alpha = 0.2)
  ## Add US states
  p <- p + geom_polygon(data = map_data("state"), aes(x = long, y = lat, group = group), 
                        colour = "black", fill = "white", alpha = 0)
  ## Modify legend and axes
  p <- p + scale_colour_discrete(name="Origin", label=c("Predicted", "True")) # labels
  p <- p + scale_x_continuous(breaks = seq(-120, -70, by = 10))
  p <- p + xlab("Longitude") + ylab("Latitude") + theme_bw()
  p <- p + theme(legend.position = "bottom", legend.text = element_text(size = 12))
  print(p)
  if (save) ggsave("all-predictions.png", path = path)
}

########## Analyze Prediction Error ##########

summarize_error <- function(err) {
  c(n = nrow(err), quantile(err[, 1], c(0.5, 0.05, 0.95)), colMeans(err[, -1]))
}

summarize_error_by_variable <- function(err, var, probs = c(0.33, 0.67)) {
  by.var <- split(err, cut(var, breaks = quantile(var, probs = unique(c(0, probs, 1)))))
  table.var <- do.call(rbind, lapply(by.var, summarize_error))
  return(table.var)
}