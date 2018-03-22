library(snpStats)

dat <- read.impute("datasim/CEU.0908.impute.files/example.gen")
snpsum.col <- col.summary(dat)

# thresholds for filtering
call <- 0.95
minor <- 0.01

use <- with(snpsum.col, (!is.na(MAF) & MAF > minor) & Call.rate >= call)
use[is.na(use)] <- FALSE

cat(ncol(dat) - sum(use), "SNP's removed due to low MAF or call rate.\n")

dat <- dat[, use]
snpsum.col <- snpsum.col[use,]

hardy = 10^-6      # HWE cut-off

snpsum.colCont = col.summary(dat)
HWEuse = with(snpsum.colCont, !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hardy/2) ) ) )
rm(snpsum.colCont)
HWEuse[is.na(HWEuse)] = FALSE          # Remove NA's as well
cat(ncol(dat)-sum(HWEuse),"SNPs will be removed due to high HWE.\n")

dat <- dat[, HWEuse]
snpsum.col <- snpsum.col[HWEuse,]

X = as(dat, "numeric")

Sigma = cov(X)
Sigma.distance = as.dist(1 - abs(cov2cor(Sigma)))
fit = hclust(Sigma.distance, method="single")
corr_max = 0.75
clusters = cutree(fit, h=1-corr_max)
