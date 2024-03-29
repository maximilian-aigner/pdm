library(knockoff)
library(glmnet)
library(SGL)
library(caret, warn.conflicts = FALSE)

stat.combined <- function(X, X_k, y, combination_function, ret.copy = FALSE) {
  out <- c()
  # out <- rbind(out, stat.random_forest(X, X_k, y))
  # out <- rbind(out, stat.stability_selection(X, X_k, y))
  # out <- rbind(out, stat.sqrt_lasso(X, X_k, y))
  out <- rbind(out, stat.glmnet_coefdiff(X, X_k, y, family = "binomial")) # lasso
  out <- rbind(out, stat.glmnet_lambdadiff(X, X_k, y, family = "binomial")) # diff of lambdas
  out <- rbind(out, stat.glmnet_coefdiff(X, X_k, y, family = "binomial", alpha = 0)) # ridge
  out <- rbind(out, stat.xgboost(X, X_k, y, n.cv = 2))
  if (ret.copy)
    return(list(Wmat = out, combined = combination_function(out)))
  return(combination_function(out))
}

stat.combined.groups <- function(X, X_k, y, groups, combination_function, ret.copy = FALSE, ...) {
  out <- c()
  out <- rbind(out, stat.glmnet_coefdiff(X, X_k, y))
  out <- rbind(out, stat.group_logit_lasso(X, X_k, y, groups, penalty = "cMCP", ...)$W)
  # out <- rbind(out, stat.group_logit_lasso(X, X_k, y, groups, penalty = "cMCP")$W)
  out <- rbind(out, stat.xgboost(X, X_k, y, n.cv = 2))
  if (ret.copy)
    return(list(Wmat = out, combined = combination_function(out)))
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
    weights <- 1.0/apply(Wmat, 1, function(row) max(row) - min(row))
  else if (weights == "var")
    weights <- 1.0/apply(Wmat, 1, stats::var)
  # don't forget to normalise weights, you dummy
  weights <- weights * 1/sum(weights)
  if (type == "sum")
    Wfinal <- colSums(weights*Wmat)
  else if (type == "mean")
    Wfinal <- colMeans(weights*Wmat)
  return(list(
    Wfinal = Wfinal,
    weights = weights,
    correlation = cor(t(Wmat))
  ));
}

highest.lambda <- function(fit, nnz) {
  lambdas <- fit$lambda
  nnz_per_lambda <- sapply(lambdas, function(l) sum(coef(fit, lambda = l) != 0))
  lambda = max(lambdas[nnz_per_lambda >= nnz])
  return(list(lambda = lambda, idx = which(lambda == lambdas)))
}

stat.group_logit_lasso <- function(X, X_k, y, groups, penalty = "grLasso", mode = "best", ...) {
  if (is.numeric(mode)) {
    if (length(mode) == 1) {
      grp.fit <- grpreg(cbind(X, X_k), y, groups, family = "binomial",
                        penalty = penalty, nlambda = 10 + mode^2, ...)
      chosen.lambda <- highest.lambda(grp.fit, mode)$lambda
      Z = coef(grp.fit, lambda = chosen.lambda)
      selection.losses = NULL
    } else {
      # mode contains several values of nnz
      # for each, compute the lowest lambda
      # select the one with least error among these lambdas
      selected.lambdas <- c()
      losses <- c()
      coef.matrix <- c()
      for (i in seq_along(mode)) {
        current.nnz <- mode[i]
        grp.fit <- grpreg(cbind(X, X_k), y, groups, family = "binomial", penalty = penalty, nlambda = 30 + current.nnz^2, ...)
        hl <- highest.lambda(grp.fit, current.nnz)
        selected.lambdas[i] <- hl$lambda
        losses[i] <- grp.fit$loss[hl$idx]
        coef.matrix <- rbind(coef.matrix, coef(grp.fit, lambda = hl$lambda))
      }
      least.loss.index <- which(losses == min(losses))
      chosen.lambda <- selected.lambdas[least.loss.index]
      Z <- coef.matrix[least.loss.index, ]
      selection.losses <- list(lambda = selected.lambdas, loss = losses)
    }
  } else {
    # assume we are choosing the best (CV sense) lambda
    grp.fit <- cv.grpreg(cbind(X, X_k), y, groups, family = "binomial",
                         penalty = penalty, ...)
    chosen.lambda <- grp.fit$lambda.min
    selection.losses <- list(lambda = grp.fit$lambda, loss = grp.fit$cve)
  }
  Z = coef(grp.fit, lambda = chosen.lambda)
  p = dim(X)[2]
  orig = 2:(p + 1)
  W = abs(Z[orig]) - abs(Z[p + orig])
  return(list(W = W, selection.losses = selection.losses))
}

stat.sparse_group_lasso <- function(X, X_k, y, groups, intercept = TRUE, mode = "best", penalty = NULL, ...) {
  data <- list(
    x = if (intercept) model.matrix(y ~ X + X_k) else model.matrix(y ~ X + X_k - 1),
    y = y
  )
  sgl.cv.fit <- cvSGL(data, type = "logit", index = groups,
               nlam = 100, standardize = FALSE, alpha = 1, ...)
  if (mode == "best") {
    least.lambda.index <- which.min(sgl.cv.fit$lldiff)
    Z <- sgl.cv.fit$fit$beta[, least.lambda.index]
  } else {
    # assume a vector of nnzs
    nnzs <- colSums(sgl.cv.fit$fit$beta != 0)
    nnz.indices <- match(mode, nnzs)
    nnz.indices <- nnz.indices[!is.na(nnz.indices)]
    least.lambda.index <- which.min(sgl.cv.fit$lldiff[nnz.indices])
    Z <- sgl.cv.fit$fit$beta[, nnz.indices[least.lambda.index]]
  }
  p <- dim(X)[2]
  orig <- 2:(p + 1)
  W <- abs(Z[orig]) - abs(Z[p + orig])
  return(list(W = W, selection.losses = sgl.cv.fit$lldiff))
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
    max_depth = c(6),
    eta = c(0.05),
    gamma = c(0.01),
    colsample_bytree = c(0.5, 0.75),
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