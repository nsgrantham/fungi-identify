
#####################################
#########  PLOT OCCURRENCE  #########
#####################################

# make sure you have run source("get-data.R") at least once!
source("set-workspace.R")  # takes some time to load Y

# Possible kernel bandwidths (defined by rho in the accompanying manuscript)
bw <- c(100, 200, 300, 400, 500, 1000, 10000)  

# Create prediction grid T across USA
# Finer resolution grid requires increased computation time
Tgrid <- generate_grid(x = 100, y = 100)

# Make occurrence maps
# Explore tax to find which fungi to plot
plot_occurrence("OTU_5851", tax, bw, Tgrid, save = TRUE)
plot_occurrence(c("OTU_232", "OTU_29971"), tax, bw, Tgrid, save = TRUE)

# Eutypa lata from the manuscript (differs slightly from Fig 1 due to new data)
# There appears to be two distinct strains of this species
plot_occurrence(rownames(tax)[which(tax$tax7 == "s__Eutypa lata")], tax, bw, Tgrid, save = TRUE)

