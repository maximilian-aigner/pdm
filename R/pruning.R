library(snpStats)

pruning.clusters <- function(X, ...) {
  clusters <- grouping.clusters(X, ...)
  ind.screen = rep(F, length(phenotypes))
  ind.screen[sample.int(length(phenotypes), size=length(phenotypes)*0.2)] = T
  pvals.screen = p.value(single.snp.tests(phenotypes[ind.screen], snp.data = genotypes[ind.screen,]), df=1)
  ind.repr = sapply(1:max(clusters), function(c) {
    cluster_elements = clusters==c
    top_within = which.min(pvals.screen[cluster_elements])
    if( length(top_within)==0 ) top_within = 1
    which(cluster_elements)[top_within]
  })
  X = X[,ind.repr]
  pvals.screen = pvals.screen[ind.repr]
  
  Xk = invisible(hmm.knockoffs(X))
  colnames(Xk) <- paste0(colnames(X), "_knockoff")
  
  Xk[ind.screen,] = X[ind.screen,]
  
  return(list(X=X, Xk=Xk, screened=ind.screen))
}