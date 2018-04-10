library(plyr)
source('./R/annotation.R')

grouping.clusters <- function(X, corr_max = 0.75, method = "single") {
  # Use clusters as groups
  Sigma <- cov(X)
  Sigma.distance <- as.dist(1 - abs(cov2cor(Sigma)))
  fit <- hclust(Sigma.distance, method=method)
  clusters <- cutree(fit, h = 1 - corr_max)
  return(clusters)
}

grouping.annotations <- function(X, legend_file="./datasim/working_dataset/hapgen2/generated_output.legend") {
  ann.df <- read.table(legend_file, sep = ' ', header = TRUE, stringsAsFactors = FALSE)
  ann.df <- cbind(chr=rep("chr22", nrow(ann.df)), ann.df[, c("rs", "pos")])
  annotations <- annotate.snps_with_genes(ann.df)
  copy <- colnames(X)
  names(copy) <- colnames(X)
  mapvalues(copy, annotations$names, annotations$GENESYMBOL)
  return(copy)
}
