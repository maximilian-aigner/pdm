grouping.clusters <- function(X, corr_max = 0.75, method = "single") {
  # Use clusters as groups
  Sigma <- cov(X)
  Sigma.distance <- as.dist(1 - abs(cov2cor(Sigma)))
  fit <- hclust(Sigma.distance, method=method)
  clusters <- cutree(fit, h = 1 - corr_max)
  return(clusters)
}

grouping.annotations <- function() {
  
}