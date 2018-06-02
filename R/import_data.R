library(snpStats, warn.conflicts = FALSE)

source('./R/utils.R')
source('./R/qc.R')

import.data <- function() {
  data.dir <- './datasim/working_dataset/hapgen2/'
  dat.cases <- read.impute(paste0(data.dir, 'generated_output.cases.gen'))
  dat.contr <- read.impute(paste0(data.dir, 'generated_output.controls.gen'))
  row.names(dat.cases) <- paste("Case", 1:dim(dat.cases)[1])
  row.names(dat.contr) <- paste("Control", 1:dim(dat.contr)[1])
  genotypes = rbind(dat.contr, dat.cases)
  phenotypes = c(rep(0,dim(dat.contr)[1]), rep(1,dim(dat.cases)[1]))
  
  idx.kept <- qc(genotypes)
  genotypes <- genotypes[, idx.kept]
  snpsum.col <- snpsum.col[idx.kept, ]
  X = as(genotypes, "numeric")
  return(list(X = X, phenotypes = phenotypes))
}