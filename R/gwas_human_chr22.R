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

set.seed(43192)

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
  groups <- grouping.annotations(X, singletons.aggregate = TRUE)
  
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
W = stat.combined.groups(X, Xk, phenotypes, total.groups, combine.weighted)
pdf("figures/combinedW_grlasso_cMCP_xgboost.pdf")
t = knockoff.threshold(W$combined, offset = 0, fdr = .15)
plot.discoveries(W$combined, t = t)
dev.off()

# W <- c()
# for (config in wanted.plots) {
#   W = rbind(W, stat.group_logit_lasso(X, Xk, phenotypes, total.groups, penalty = config$penalty, mode = config$mode))
# }

#pdf("figures/combinedW_aggregated.pdf")
#  par(cex.main = 1.2, cex.lab = 1.2)
#  par(mar = c(5.1, 5.1, 3.1, 2.1))
#  ncases <- ceiling(length(wanted.plots)/2)*2
#  layout(matrix(1:ncases, ncol = 2, byrow = T))
#  for (i in 1:length(wanted.plots)) {
#    t = knockoff.threshold(W[i, ]$W, offset = 0, fdr = .15)
#    plot.discoveries(W[i, ]$W, t = t,
                   #main = penalty.names[[wanted.plots[[i]]$penalty]], 
#                   ylim = c(-.3, .2))
#  }
#dev.off()

# Wmat.combined <- stat.combined.groups(X, Xk, phenotypes, total.groups, combine.weighted, ret.copy = TRUE)
# W.combined <- Wmat.combined$combined

# Compare to true active set
active <- read.table('~/src/pdm/datasim/working_dataset/active_genes.txt',
                     stringsAsFactors = FALSE, header = TRUE, sep = ' ')
active.genes <- unique(active$GENESYMBOL)
active.genes.snps <- which(groups %in% active.genes)
names(active.genes.snps) <- names(groups[active.genes.snps])