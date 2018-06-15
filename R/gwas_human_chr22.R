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
  grouping.obj <- grouping.clusters(X)
  groups <- grouping.obj$clusters
  X <- grouping.obj$Xmod
  
  save(X, groups, file = './preload/savegroups.rda')
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

wanted.plots <- list(
  list(
    f = stat.group_logit_lasso,
    penalty = "cMCP",
    mode = "best"
  ),
  list(
    f = stat.group_logit_lasso,
    penalty = "cMCP",
    mode = 1:20
  )
)

W <- c()
for (config in wanted.plots) {
  cat("[plotting] (", config$penalty, ", ", config$mode, ")\n")
  output <- do.call(
    config$f, 
    list(X = X, X_k = Xk,
         y = phenotypes, groups = total.groups,
         penalty = config$penalty, mode = config$mode
         )
  )
  W = rbind(W, output$W)
}

pdf("figures/cluster_bilevelselection.pdf")
 par(cex.main = 1.2, cex.lab = 1.2)
 par(mar = c(5.1, 5.1, 3.1, 2.1))
 ncases <- ceiling(length(wanted.plots)/2)*2
 layout(matrix(1:ncases, ncol = 1, byrow = T))
 for (i in 1:length(wanted.plots)) {
   t = knockoff.threshold(W[i, ], offset = 0, fdr = .15)
   plot.discoveries(W[i, ], t = t,
                  ylim = c(-1, 2))
 }
dev.off()

# Compare to true active set
active <- read.table('~/src/pdm/datasim/working_dataset/active_genes.txt',
                     stringsAsFactors = FALSE, header = TRUE, sep = ' ')
active.genes <- unique(active$GENESYMBOL)
active.genes.snps <- which(groups %in% active.genes)
names(active.genes.snps) <- names(groups[active.genes.snps])

# indirect.good <- intersect(names(disc.obj$discoveries), names(active.genes.snps))
# actual.fdp <- length(indirect.good) / length(disc.obj$discoveries)
