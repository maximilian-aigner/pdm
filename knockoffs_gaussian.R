source('utils.R')
library(pracma)

orthog.basis <- function(X) {
  
}

project <- function(u, v) {
  up <- (u %*% v)/(v %*% v);
  return(as.numeric(up))
}

gramschmidt <- function(df, indices) {
  # Orthogonalise those columns corresponding to indices
  sz <- dim(df)
  p <- sz[2]
  fprintf("Orthogonalising columns:\n")
  print(indices)
  fprintf("Among 1:%d total\n", p)
  alli <- 1:p
  is_comp <- alli[!(alli %in% indices)]
  for (i in indices)
  {
    fprintf("Doing column %d\n", i)
    # orthogonalise against all original cols,
    # plus those that have already been entered
    for (j in is_comp)
    {
      fprintf("Against column %d\n", j)
      prj <- project(df[, i], df[, j])
      fprintf("Projection = %.5f\n", prj)
      df[, i] <- df[, i] - prj * df[, j]
    }
    #df[, i] <- df[, i] - sum(
    #  apply(df[, is_comp], 2, function(u) project(u, df[, i]))
    #)
    is_comp <- c(is_comp, i)
  }
  return(df)
}

knockoffs.g <- function(X) {
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
  Xtilde <- X %*% (I - sig.inv) - U %*% C
  aug.X <- cbind(X, Xtilde)
  return(aug.X)
}

