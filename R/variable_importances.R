library(knockoff)
library(caret)

stat.combined <- function(X, X_k, y, combination_function) {
  out <- c()
  #out <- rbind(out, stat.random_forest(X, X_k, y))
  #out <- rbind(out, stat.stability_selection(X, X_k, y))
  out <- rbind(out, stat.sqrt_lasso(X, X_k, y))
  out <- rbind(out, stat.glmnet_coefdiff(X, X_k, y))
  out <- rbind(out, stat.glmnet_lambdadiff(X, X_k, y))
  out <- rbind(out, stat.glmnet_coefdiff(X, X_k, y, alpha = 0))
  out <- rbind(out, stat.xgboost(X, X_k, y))
  return(combination_function(out))
}

combine.max <- function(Wmat) {
  W <- apply(Wmat, 2, max)
  return(W)
}

combine.prod <- function(Wmat) {
  W <- apply(Wmat, 2, prod)
}

combine.weighted <- function(Wmat, type = "sum", weights = "sd") {
  if (weights == "sd")
    weights <- 1.0/apply(Wmat, 1, stats::sd)
  else if (weights == "range")
    weights <- 1.0/apply(Wmat, 1, function(row) max(row)-min(row))
  if (type == "sum")
    Wfinal <- colSums(Wmat*weights)
  else if (type == "mean")
    Wfinal <- colMeans(Wmat*weights)
  return(Wfinal);
}

stats.group_logit_lasso <- function(X, X_k, y, groups, penalty = "grLasso", mode = "best", ...) {
  if (is.numeric(mode)) {
    # minimum guaranteed nonzero coefs (mode = number of them)
    # heuristic: number of lambdas should be at least mode^2
    min.vals <- mode
    min.vals.lambdas <- c()
    min.vals.errs <- c()
    for (i in 1:length(min.vals)) {
      min_coefs <- min.vals[i]
      one.grp.fit <- grpreg(cbind(X, X_k), y, groups, family = "binomial",
                         penalty=penalty, nlambda = 1+min_coefs^2, ...)
      grp.lambdas <- one.grp.fit$lambda
      nzcoefs <- sapply(grp.lambdas, function(l) sum(coef(one.grp.fit, lambda = l) != 0))
      highest.lambda <- max(grp.lambdas[nzcoefs >= min_coefs])
      highest.lambda.idx <- which(grp.lambdas[nzcoefs >= min_coefs] == highest.lambda)
      min.vals.lambdas[i] <- highest.lambda
      min.vals.errs[i] <- one.grp.fit$loss[highest.lambda.idx]
    }
    # now, select the (lambda, nnz) with lowest loss
    best.overall <- which(min.vals.errs == min(min.vals.errs))
    chosen.lambda <- min.vals.lambdas[best.overall]
    print(min.vals.lambdas)
    grp.fit <- grpreg(cbind(X, X_k), y, groups, family = "binomial", penalty = penalty, lambda = min.vals.lambdas, max.iter = 1e6)
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

stat.xgboost <- function(X, X_k, y, n.cv = 4) {
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

