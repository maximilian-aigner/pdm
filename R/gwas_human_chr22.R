library(snpStats)
library(knockoff)
library(grpreg)
library(glmnet)
library(plyr)

source('./R/variable_importances.R')
source('./R/utils.R')
source('./R/grouping.R')
source('./R/qc.R')

set.seed(43192)

# Read in files
dat.cases <- read.impute("./datasim/working_dataset/hapgen2/generated_output.cases.gen")
dat.contr <- read.impute("./datasim/working_dataset/hapgen2/generated_output.controls.gen")
row.names(dat.cases) = sapply(1:dim(dat.cases)[1], function(i) paste("Case",i,sep=""))
row.names(dat.contr) = sapply(1:dim(dat.contr)[1], function(i) paste("Control",i,sep=""))
genotypes = rbind(dat.contr, dat.cases)
phenotypes = c(rep(0,dim(dat.contr)[1]), rep(1,dim(dat.cases)[1]))
snpsum.col <- col.summary(genotypes)

# Quality controls
idx.kept <- qc(genotypes)
genotypes <- genotypes[, idx.kept]
snpsum.col <- snpsum.col[idx.kept,]

# Now, analyse the remaining columns
X = as(genotypes, "numeric")

# Generate groups by clustering
groups <- grouping.annotations(X)

# Generate knockoffs
Xk = invisible(hmm.knockoffs(X))

# Group knockoffs as the originals, but not paired with them
total.groups <- c(groups, paste0(groups, "_knockoff"))
names(total.groups) <- c(names(groups), paste0(names(groups), "_knockoff"))


W = stats.group_logit_lasso(X, Xk, phenotypes, total.groups)
thresh = knockoff.threshold(W, fdr = 0.1, offset = 0)

plot.discoveries(W, thresh)

