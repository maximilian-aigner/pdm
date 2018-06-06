library(SNPknock, warn.conflicts = FALSE)
library(scales)

fprintf <- function(fmt, ..., file = "", append = FALSE) {
  mystr <- sprintf(fmt, ...)
  cat(mystr, file = file, append = append)
  invisible(nchar(mystr))
}

hmm.knockoffs <- function(X, ...) {
  # run fastPHASE
  # Xout_path <- "./datasim/fastPHASE/OUTPUT_fastPHASE/X_input"
  # fp_output_path = "./datasim/fastPHASE/OUTPUT_fastPHASE"
  # fp_outPath = SNPknock.fp.runFastPhase(fp_path, Xinp_file, out_path = fp_output_path, K = 8, numit = 10, ...)
  
  storage.mode(X) <- "integer"
  Xinp_file = SNPknock.fp.writeX(X)
  fp_path  = "./datasim/fastPHASE/fastPHASE.bin"
  fp_outPath = SNPknock.fp.runFastPhase(fp_path, Xinp_file, K = 8, numit = 10)
  
  r_file = paste(fp_outPath, "_rhat.txt", sep="")
  theta_file = paste(fp_outPath, "_thetahat.txt", sep="")
  alpha_file = paste(fp_outPath, "_alphahat.txt", sep="")
  char_file = paste(fp_outPath, "_origchars", sep="")
  hmm = SNPknock.fp.loadFit(r_file, theta_file, alpha_file, char_file)
  
  Xk = SNPknock.knockoffGenotypes(X, hmm$r, hmm$alpha, hmm$theta)
  return(Xk)
}

plot.discoveries <- function(W, t, plotit = TRUE, ...) {
  discoveries = which(W >= t)
  names(discoveries) = names(W)[discoveries]

  real <- get.active.snps()
  signals <- real$rs
  signals.id <- match(signals, names(W))
  names(signals.id) <- signals
  
  if (plotit) { 
    plot(W, col = alpha("gray", .5), pch = 16, cex = 1, ylab = expression('W'[j]), xlab = 'j', ...) 
    abline(h = t, lty = 2)
    points(signals.id, W[signals.id], col = "#377EB8", pch = 0)
    points(discoveries, W[discoveries], col = "#E41A1C", pch = 4)
  }
  return(list(discoveries = discoveries, signals = signals.id))
}

empirical.fdr <- function(W, t, ...) {
  outcomes <- plot.discoveries(W, t, plotit = FALSE, ...)
  false <- setdiff(outcomes$signals, outcomes$discoveries)
  empirical <- length(false)/length(discoveries)
  return(empirical)
}

get.active.snps <- function(fname = "datasim/working_dataset/active_genes.txt") {
  real <- read.table(fname, header = TRUE, stringsAsFactors = FALSE)
  return(real)
}
