source('utils.R')
library(pracma)

knockoffs.g <- function(X) {
  mu <- colMeans(X)
  sig <- cov(X)
  sz <- dim(X)
  n <- sz[1]
  p <- sz[2]

  # equicorrelated
  lambda_0 <- eigen(sig)$values[p]

  s <- 2 * rep(1, p) * lambda_0
  S <- diag(s)
  sig.inv <- solve(sig, S)

  C <- chol(2 * S - S * sig.inv)
  U <- matrix(rnorm(n*p), ncol = p)

  I = diag(1, p)
  # mut <- X - sweep(X, 2, mu, "-") %*% sig.inv
  
  Xt <- X %*% (I - sig.inv) - U %*% C
  #Xt <- mut + U %*% C
  aug.X <- cbind(X, Xt)
  return(aug.X)
}

