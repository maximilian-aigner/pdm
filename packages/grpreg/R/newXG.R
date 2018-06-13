newXG <- function(X, g, m, ncolY, bilevel) {
  # Coerce X to matrix
  if (class(X) != "matrix") {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (any(is.na(X))) stop("Missing data (NA's) detected in X.  You must eliminate missing data (e.g., by removing cases, removing features, or imputation) before passing X to grpreg")

  # Reorder groups, if necessary
  if (length(g) != ncol(X)) stop ("Dimensions of group is not compatible with X")
  xnames <- if (is.null(colnames(X))) paste("V",1:ncol(X),sep="") else colnames(X)
  grp <- reorderGroups(g, m, bilevel)
  g <- grp$g
  m <- grp$m
  if (grp$reorder) X <- X[,grp$ord]

  # Make multiX, if necessary
  if (ncolY > 1) {
    X <- multiX(X, ncolY)
    g <- c(rep(0, ncolY-1), rep(g, each=ncolY))
  }

  # Standardize
  std <- .Call("standardize", X)
  XX <- std[[1]]
  cat("[newXG.R] dimensions of X:", dim(X), "\n")
  cat("[newXG.R] dimensions of XX:", dim(XX), "\n")
  center <- std[[2]]
  scale <- std[[3]]
  if (anyNA(scale))
  {
    the.nas <- which(is.na(scale))
    cat("[newXG.R] replacing %d variables by N(0, 1)\n", length(the.nas))
    XX[, the.nas] <- rnorm(nrow(XX))
    scale[the.nas] <- 1
  }
  nz <- which(scale > 1e-6)                # non-constant columns
  if (any(is.na(nz)))
    stop('NA SCALE')
  zg <- setdiff(unique(g), unique(g[nz]))  # constant groups
  if (length(zg)) {
    cat("[newXG.R] Found constant groups (length(zg) = ", length(zg), ") (will be imputed)\n")
    gf <- factor(g)
    #gf <- factor(g[!(g %in% zg)])
    if (any(levels(gf)=="0")) {
      cat("[newXG.R] found 0-encoded factor level\n")
      g <- as.numeric(gf) - 1
    } else {
      g <- as.numeric(gf)
    }
    # m <- m[-zg]
  }
  if (length(nz) != ncol(X)) {
    cat("[newXG.R] imputing ")
    problems<-setdiff(1:ncol(X), nz)
    cat(length(problems), " constant variable(s)\n")
    # print(problems)
    # cat("because scales=")
    # print(scale[problems])
    scale[problems] <- 1
    XX[, problems] <- rnorm(length(problems) * nrow(XX))
    nz <- 1:ncol(X) 
    # XX <- XX[ ,nz, drop=FALSE]
    # change g also
    # cat("`nz` NAs: ")
    # print(nz[is.na(nz)])
    # save(nz, g, file = '~/src/pdm/debug/nz_and_g.Rda')
    # cat("before removal, g has", sum(is.na(g)), "NAs\n")
    # g <- g[nz]
    # cat("after removal, g has", sum(is.na(g)), "NAs\n")
  }
  if (!bilevel) {
    XX <- orthogonalize(XX, g)
    g <- attr(XX, "group")
  }

  # Return
  return(list(X=XX, g=g, m=m, reorder=grp$reorder, ord.inv=grp$ord.inv, names=xnames,
              center=center[nz], scale=scale[nz], nz=nz))
}
