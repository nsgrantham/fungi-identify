
#####################################
############  GET DATA  #############
#####################################

# Data files stored at
# http://figshare.com/articles/1000homes/1270900/

dir.create("data")
dir.create("data/raw")

## download homes txt file (200 Kb)
download.file("https://ndownloader.figshare.com/files/3254927", "data/raw/homes_mapping_file.txt")

## download biom files (~70 Mb)
# create a temporary directory
td <- tempdir()
# create the placeholder file
tf <- tempfile(tmpdir=td, fileext=".zip")
# download into the placeholder file
download.file("https://ndownloader.figshare.com/files/3254933", tf)
# unzip the file to the raw directory
unzip(tf, files = "ITS_otu_table_wTax.biom", exdir = "data/raw", overwrite = TRUE)
# delete temp file
unlink(tf)

## download fa files (~5 Mb) 
# (NOTE: These files are not necessary for the statistical analyses)
tf <- tempfile(tmpdir=td, fileext=".zip")
download.file("https://ndownloader.figshare.com/files/3254954", tf)
unzip(tf, files = "ITS_rep_set_numbered.fa", exdir = "data/raw", overwrite = TRUE)
unlink(c(tf, td))
rm(list = c("tf", "td"))

##################################
##########  MUNGE DATA  ##########
##################################

# Brief description of the following data:
# X -- covariate data related to home dust samples
# S -- lon, lat of sample coordinates (NOTE: exact coordinates have been 
#      truncated to protect confidentiality of study participants.)
# Y -- presence/absence of taxa by DNA sequencing of dust samples
# tax -- taxonomic designation of fungal taxa

#####  X  #####

X <- read.table("data/raw/homes_mapping_file.txt", header=TRUE, sep="\t")
X <- X[!duplicated(X[, 1]), ]  # remove one or two duplicate entries
# make home IDs rownames of X
rownames(X) <- X[, 1]
X <- X[, -1]
# since authoring of the manuscript, some homes have been added from AK and HI.
# To preserve continuity with the manuscript, we remove them
X <- X[!(X[, "State"] %in% c("AK", "HI")), ]

#####  S  #####
# will be helpful to separate lon lat coordinates from X

S <- X[, c("Longitude", "Latitude")]
colnames(S) <- c("lon", "lat")
X <- X[, -which(colnames(X) %in% c("Longitude", "Latitude"))]

#####  Y  #####

# The biom package has not been updated in several years and was removed from CRAN.
# Therefore, `install.packages("biom")` will not work. We must install an archived version.
# biom depends on RJSONIO, so install it before we begin with `install.packages("RJSONIO")`
# Next Go to https://cran.r-project.org/src/contrib/Archive/biom/ and download the biom_0.3.12.tar.gz.
# Navigate to the directory containing biom_0.3.12.tar.gz in your terminal and run:
# RMD CMD INSTALL biom_0.3.12.tar.gz.
library(biom)
fungi.biom <- read_biom("data/raw/ITS_otu_table_wTax.biom")

# WARNING: following object is large (~1.5 Gb)
Y <- t(as.data.frame(as(biom_data(fungi.biom), "matrix")))
# rownames of Y are each home ID and sample type (I = Indoor, O = Outdoor)
IDtype <- strsplit(rownames(Y), "\\.")
# retain only outdoor samples for each home
is.outdoor <- sapply(IDtype, function(x) x[2] == "O")
# also remove samples that did not come from homes (ExtB1, NTC12, ...)
do.call(c, IDtype[is.na(is.outdoor)])
is.home <- !is.na(is.outdoor)
Y <- Y[is.home & is.outdoor, ]
IDtype <- IDtype[is.home & is.outdoor]
# now order home IDs to match rows of X & S
ID <- as.numeric(sapply(IDtype, function(x) x[1]))  # some NAs introduced
Y <- Y[!is.na(ID), ]
ID <- ID[!is.na(ID)]
Y <- Y[order(ID), ]
rownames(Y) <- ID[order(ID)]

# remove homes from X & S that do not appear in Y
# This occurs when an outdoor dust sample failed to sequence correctly
X <- X[rownames(X) %in% rownames(Y), ]
S <- S[rownames(S) %in% rownames(Y), ]
all(rownames(X) == rownames(S))  # S and X have same homes
# there is one home in Y that is not in X or S (why?), so just remove it
Y <- Y[rownames(Y) %in% rownames(X), ]
# do all IDs match now?
all(rownames(X) == rownames(Y))  # success!

# currently Y is number of reads for each OTU (col) & each home (row)
# for simplicity, change Y to binary. i.e., taxa is present vs. absent
Y[Y > 0] <- 1
# Some columns in Y are entirely 0, so taxa was not sequenced in any outdoor dust sample.
# (However, the taxa was sequenced from another sample type in the study)
# Remove these fungi taxa that are never observed in the data
use <- colSums(Y) > 0
Y <- Y[, use]

# Quick look at presence/absence fungal taxa data
dim(Y)  # n-by-m
mean(Y)  # very sparse
hist(rowSums(Y), breaks = 20)  # number of taxa sequenced per sample
mean(rowSums(Y))  # avg. number of taxa per sample
mean(colSums(Y) < 10)   # proportion of taxa found in <10 of n samples
mean(colSums(Y) < 100)  # proportion of taxa found in <100 of n samples

#####  tax  #####

fungi.tax <- observation_metadata(fungi.biom)
tax <- data.frame(tax1 = sapply(fungi.tax, function(x) x[1]),
                  tax2 = sapply(fungi.tax, function(x) x[2]),
                  tax3 = sapply(fungi.tax, function(x) x[3]),
                  tax4 = sapply(fungi.tax, function(x) x[4]),
                  tax5 = sapply(fungi.tax, function(x) x[5]),
                  tax6 = sapply(fungi.tax, function(x) x[6]),
                  tax7 = sapply(fungi.tax, function(x) x[7]))
rownames(tax) <- names(fungi.tax)
tax <- tax[use, ]  # restrict to only those taxa observed in Y

# do the OTUs in tax match those in Y?
all(rownames(tax) == colnames(Y))  # success!

# finally, write to csv
write.csv(S, "data/S.csv")
write.csv(X, "data/X.csv")
write.csv(Y, "data/Y.csv")
write.csv(tax, "data/tax.csv")

# in the interest of file size, share a versionof the data with
# only the taxa that appear more than 100 samples
is.prevalent <- colSums(Y) > 100
write.csv(Y[, is.prevalent], "data/Y-reduced.csv")
write.csv(tax[is.prevalent, ], "data/tax-reduced.csv")

# and remove variables from workspace
rm(list = ls())