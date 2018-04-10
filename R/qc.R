# qc.R: perform standard quality controls
qc <- function(genotypes, call = 0.95, maf.minor = 0.01, hwe = 1e-6, verbose = FALSE) {
  use <- with(snpsum.col, (!is.na(MAF) & MAF > maf.minor) & Call.rate >= call)
  use[is.na(use)] <- FALSE
  
  if(verbose)
    cat(ncol(genotypes) - sum(use), "SNP's removed due to low MAF or call rate.\n")
  
  genotypes <- genotypes[, use]
  
  snpsum.colCont = col.summary(genotypes)
  HWEuse = with(snpsum.colCont, !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hwe/2) ) ) )
  rm(snpsum.colCont)
  
  HWEuse[is.na(HWEuse)] = FALSE
  if(verbose)
    cat(ncol(genotypes)-sum(HWEuse), "SNP's will be removed due to high HWE.\n")
  padded <- rep(FALSE, length(use))
  padded[HWEuse] <- TRUE
  final.idx <- use & padded
  return(final.idx)
}