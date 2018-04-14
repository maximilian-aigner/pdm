library(plyr)
source('./R/annotation.R')

grouping.clusters <- function(X, corr_max = 0.75, method = "single", draw.mode = 0, ...) {
  # Use clusters as groups
  Sigma <- cov(X)
  cov.cor <- cov2cor(Sigma)
  Sigma.distance <- as.dist(1 - abs(cov.cor))
  fit <- hclust(Sigma.distance, method = method)
  clusters <- cutree(fit, h = 1 - corr_max)
  # Appropriate plots
  if (draw.mode == 1) {
    library("heatmap3")
    heatmap3(cov.cor, useRaster = TRUE, ...)
  } else if (draw.mode == 2) {
    hist(log10(clusters), xlim = rev(range(log10(clusters))),
         xlab = "log(group size)", ylab = "", col = "gray", main = "Cluster sizes", ...)
  }
  return(clusters)
}

grouping.annotations <- function(X, legend_file="./datasim/working_dataset/hapgen2/generated_output.legend", verbose = FALSE) {
  ann.df <- read.table(legend_file, sep = ' ', header = TRUE, stringsAsFactors = FALSE)
  ann.df <- cbind(chr = rep("chr22", nrow(ann.df)), ann.df[, c("rs", "pos")])
  annotations <- annotate.snps_with_genes(ann.df)
  annotations <- annotations[complete.cases(annotations), ]
  copy <- colnames(X)
  names(copy) <- colnames(X)
  copy <- mapvalues(copy, annotations$names, annotations$GENESYMBOL, warn_missing = FALSE)
  if (verbose) {
    n.unchanged <- sum(names(copy) == copy)
    n.changed <- length(copy) - n.unchanged
    rat.changed <- n.changed/length(copy)
    cat("Annotated", n.changed, "SNPs — ", rat.changed*100, "%.\n")
  }
  return(copy)
}
