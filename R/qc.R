# qc.R: perform standard quality controls

source('./R/utils.R')

qc <- function(genotypes, call = 0.95, maf.minor = 0.01, hwe = 1e-6, verbose = FALSE) {
  active <- get.active.snps()
  active.names <- active$rs
  
  orig.names <- colnames(genotypes)
  
  snpsum.col <- col.summary(genotypes)
  use <- with(snpsum.col, (!is.na(MAF) & MAF > maf.minor) & Call.rate >= call)
  use[is.na(use)] <- FALSE
  
  if (verbose)
    cat(ncol(genotypes) - sum(use), "SNP's removed due to low MAF or call rate.\n")
  
  genotypes <- genotypes[, use]
  
  snpsum.colCont = col.summary(genotypes)
  HWEuse = with(snpsum.colCont, !is.na(z.HWE) & (abs(z.HWE) < abs(qnorm(hwe/2))))
  rm(snpsum.colCont)
  
  HWEuse[is.na(HWEuse)] = FALSE
  if (verbose)
    cat(ncol(genotypes) - sum(HWEuse), "SNP's will be removed due to high HWE.\n")

  genotypes <- genotypes[, HWEuse]
  
  matched <- match(colnames(genotypes), orig.names)
  matched.active <- match(active.names, orig.names)
  final.idx <- rep(FALSE, length(orig.names))
  final.idx[matched] <- TRUE
  final.idx[matched.active] <- TRUE
  return(final.idx)
}