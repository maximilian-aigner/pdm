library(knockoff)
run_VIs <- function(X, X_k, y, combination_function)
{
  out <- c()
  #out <- rbind(out, stat.random_forest(X, X_k, y))
  #out <- rbind(out, stat.stability_selection(X, X_k, y))
  out <- rbind(out, stat.sqrt_lasso(X, X_k, y))
  out <- rbind(out, stat.glmnet_coefdiff(X, X_k, y))
  out <- rbind(out, stat.glmnet_lambdadiff(X, X_k, y))
  out <- rbind(out, stat.glmnet_coefdiff(X, X_k, y, alpha = 0))
  return(combination_function(out))
}

max_combination <- function(Wmat)
{
  W <- apply(Wmat, 2, max)
  return(W)
}

prod_combination <- function(Wmat)
{
  W <- apply(Wmat, 2, prod)
}

