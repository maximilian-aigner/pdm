library(knockoff)
combine_VIs <- function(X, X_k, y, combination_function)
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

stats.grouped_logit_lasso <- function(X, X_k, y, groups, penalty = "grLasso", mode = "best", ...)
{
  # grp.fit <- glmnet(cbind(X, Xk), phenotypes, family = "binomial", nlambda = 500)
  if (is.numeric(mode)) {
    # minimum guaranteed nonzero coefs (mode = number of them)
    # heuristic: number of lambdas should be at least 2*mode
    grp.fit <- grpreg(cbind(X, X_k), y, groups, family = "binomial",
                         penalty=penalty, nlambda = 2*mode, ...)
    grp.lambdas <- grp.fit$lambda
    min_coefs <- mode
    n_nzcoefs <- sapply(grp.lambdas, function(l) sum(coef(grp.fit, lambda = l)!=0))
    chosen.lambda <- max(grp.lambdas[n_nzcoefs>=min_coefs])
  } 
  else {
    # assume we are choosing the best (CV sense) lambda
    grp.fit <- cv.grpreg(cbind(X, X_k), y, groups, family = "binomial",
                         penalty=penalty, ...)
    chosen.lambda <- grp.fit$lambda.min
  }
  Z = coef(grp.fit, lambda = chosen.lambda)
  p = dim(X)[2]
  orig = 2:(p+1)
  W = abs(Z[orig]) - abs(Z[p+orig])
  return(W)
}

