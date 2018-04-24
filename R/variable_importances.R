library(knockoff)
library(caret)

combine_VIs <- function(X, X_k, y, combination_function) {
  out <- c()
  #out <- rbind(out, stat.random_forest(X, X_k, y))
  #out <- rbind(out, stat.stability_selection(X, X_k, y))
  out <- rbind(out, stat.sqrt_lasso(X, X_k, y))
  out <- rbind(out, stat.glmnet_coefdiff(X, X_k, y))
  out <- rbind(out, stat.glmnet_lambdadiff(X, X_k, y))
  out <- rbind(out, stat.glmnet_coefdiff(X, X_k, y, alpha = 0))
  return(combination_function(out))
}

max_combination <- function(Wmat) {
  W <- apply(Wmat, 2, max)
  return(W)
}

prod_combination <- function(Wmat) {
  W <- apply(Wmat, 2, prod)
}

stats.group_logit_lasso <- function(X, X_k, y, groups, penalty = "grLasso", mode = "best", ...) {
  if (is.numeric(mode)) {
    # minimum guaranteed nonzero coefs (mode = number of them)
    # heuristic: number of lambdas should be at least 5*mode
    grp.fit <- grpreg(cbind(X, X_k), y, groups, family = "binomial",
                         penalty=penalty, nlambda = mode^2, ...)
    grp.lambdas <- grp.fit$lambda
    min_coefs <- mode
    n_nzcoefs <- sapply(grp.lambdas, function(l) sum(coef(grp.fit, lambda = l)!=0))
    chosen.lambda <- max(grp.lambdas[n_nzcoefs>=min_coefs])
  } else {
    # assume we are choosing the best (CV sense) lambda
    grp.fit <- cv.grpreg(cbind(X, X_k), y, groups, family = "binomial",
                         penalty=penalty, ...)
    chosen.lambda <- grp.fit$lambda.min
  }
  Z = coef(grp.fit, lambda = chosen.lambda)
  p = dim(X)[2]
  orig = 2:(p + 1)
  W = abs(Z[orig]) - abs(Z[p + orig])
  return(W)
}

stats.xgboost <- function(X, X_k, y, n.cv = 4) {
  xgb.control <- trainControl(
    method = 'cv',
    number = n.cv,
    summaryFunction = twoClassSummary,
    classProbs = TRUE,
    verboseIter = TRUE,
    allowParallel = TRUE
  )
  xgb.grid <- expand.grid(
    nrounds = c(350),
    max_depth = c(4, 6),
    eta = c(0.05, 0.1),
    gamma = c(0.01),
    colsample_bytree = c(0.5, 0.75, 0.9),
    subsample = c(0.5),
    min_child_weight = c(0)
  )
  xgb.fit <- train(x = cbind(X, X_k), y = factor(y, labels = c("control", "case")), method = 'xgbTree', metric = 'ROC',
                   trControl = xgb.control, tuneGrid = xgb.grid)
  xgb.vi <- varImp(xgb.fit)$importance
  Zt = as.vector(t(xgb.vi))
  names(Zt) <- rownames(xgb.vi)
  Z <- c(Zt[colnames(X)], Zt[colnames(X_k)])
  p = dim(X)[2]
  orig = 1:p
  W = abs(Z[orig]) - abs(Z[p + orig])
  return(W)
}

