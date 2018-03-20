knockoff.threshold <- function(W, alpha) {
  aW <- abs(W)
  passes_threshold <- function(t) ((1 + sum(W < -t)) / sum(W > t)  < alpha)
  # print(aW[sapply(aW, passes_threshold)])
  min(aW[sapply(aW, passes_threshold)])
}