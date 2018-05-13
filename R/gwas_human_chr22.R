library(snpStats)
library(knockoff)
library(grpreg)
library(glmnet)
library(plyr)

source('./R/variable_importances.R')
source('./R/utils.R')
source('./R/grouping.R')
source('./R/qc.R')

# set.seed(43192)

# Read in files
data.dir <- './datasim/working_dataset/hapgen2/'
dat.cases <- read.impute(paste0(data.dir, 'generated_output.cases.gen'))
dat.contr <- read.impute(paste0(data.dir, 'generated_output.controls.gen'))
row.names(dat.cases) <- paste("Case", 1:dim(dat.cases)[1])
row.names(dat.contr) <- paste("Control", 1:dim(dat.contr)[1])
genotypes = rbind(dat.contr, dat.cases)
phenotypes = c(rep(0,dim(dat.contr)[1]), rep(1,dim(dat.cases)[1]))
snpsum.col <- col.summary(genotypes)

# Quality controls
idx.kept <- qc(genotypes)
genotypes <- genotypes[, idx.kept]
snpsum.col <- snpsum.col[idx.kept, ]


# Now, analyse the remaining columns
X = as(genotypes, "numeric")

# Generate groups by clustering
groups <- grouping.annotations(X)

# Generate knockoffs
Xk = invisible(hmm.knockoffs(X))
colnames(Xk) <- paste0(colnames(X), "_knockoff")

# Group knockoffs as the originals, but not paired with them
total.groups <- c(groups, paste0(groups, "_knockoff"))
names(total.groups) <- c(names(groups), paste0(names(groups), "_knockoff"))
total.groups <- as.factor(total.groups)

# Compute W-statistic
# W = stat.group_logit_lasso(X, Xk, phenotypes, total.groups, penalty = "cMCP", mode = 20)
Wo <- stat.combined.groups(X, Xk, phenotypes, total.groups, combine.prod, ret.copy = TRUE)
W <- Wo$combined

# Threshold
thresh <- knockoff.threshold(W, fdr = 0.1, offset = 0)

names(W) <- colnames(X)
outcomes <- plot.discoveries(W, thresh)

# Compare to true active set
active <- read.table('~/src/pdm/datasim/working_dataset/active_genes.txt',
                     stringsAsFactors = FALSE, header = TRUE, sep = ' ')
active.genes <- unique(active$GENESYMBOL)
active.genes.snps <- which(groups %in% active.genes)
names(active.genes.snps) <- names(groups[active.genes.snps])