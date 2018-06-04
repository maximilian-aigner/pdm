library(snpStats, warn.conflicts = FALSE)
library(knockoff, warn.conflicts = FALSE)
library(grpreg, warn.conflicts = FALSE)
library(glmnet, warn.conflicts = FALSE)
library(plyr, warn.conflicts = FALSE)

source('./R/variable_importances.R')
source('./R/utils.R')
source('./R/grouping.R')
source('./R/qc.R')
source('./R/import_data.R')

# set.seed(43192)

if (!file.exists('./preload/saveX.rda')) {
  # Re-do data import step
  loaded <- import.data()
  X <- loaded$X
  phenotypes <- loaded$phenotypes
  
  # Save the R objects for future
  save(X, phenotypes, file = './preload/saveX.rda')
} else {
  load('./preload/saveX.rda')
}

if (!file.exists('./preload/savegroups.rda')) {
  # Generate groups by clustering
  groups <- grouping.annotations(X, singletons.aggregate = FALSE)
  
  save(groups, file = './preload/savegroups.rda')
} else {
  load('./preload/savegroups.rda')
}

if (!file.exists('./preload/saveXk.rda')) {
  # Generate knockoffs
  Xk = invisible(hmm.knockoffs(X))
  colnames(Xk) <- paste0(colnames(X), "_knockoff")
  
  # Save for future reference
  save(Xk, file = './preload/saveXk.rda')
} else {
  load('./preload/saveXk.rda')
}

# Group knockoffs as the originals, but not paired with them
total.groups <- c(groups, paste0(groups, "_knockoff"))
names(total.groups) <- c(names(groups), paste0(names(groups), "_knockoff"))
total.groups <- as.factor(total.groups)

# Compute W-statistic
wanted.plots <- list(
  list(
    penalty = "grLasso", mode = "best"  
  ),
  list(
    penalty = "grLasso", mode = 1:20
  ),
  list(
    penalty = "cMCP", mode = "best"
  ),
  list(
    penalty = "cMCP", mode = 10:20
  )
)

W <- c()
for (config in wanted.plots) {
  W = rbind(W, stat.group_logit_lasso(X, Xk, phenotypes, total.groups, penalty = config$penalty, mode = config$mode))
}

layout(matrix(1:4, nrow = 2, byrow = T))
for (i in 1:3) {
  t = knockoff.threshold(W[i, ], offset = 0, fdr = .1)
  plot.discoveries(W[i, ], t = t, main = wanted.plots[i]$penalty)
}
# Wmat.combined <- stat.combined.groups(X, Xk, phenotypes, total.groups, combine.weighted, ret.copy = TRUE)
# W.combined <- Wmat.combined$combined

# Threshold
thresh <- knockoff.threshold(W, fdr = 0.5, offset = 0)

names(W) <- colnames(X)
outcomes <- plot.discoveries(W, thresh)

# Compare to true active set
active <- read.table('~/src/pdm/datasim/working_dataset/active_genes.txt',
                     stringsAsFactors = FALSE, header = TRUE, sep = ' ')
active.genes <- unique(active$GENESYMBOL)
active.genes.snps <- which(groups %in% active.genes)
names(active.genes.snps) <- names(groups[active.genes.snps])