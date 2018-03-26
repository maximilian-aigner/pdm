library(snpStats)

dat.cases <- read.impute("datasim/OUTPUT_IMPUTE2/imputed.cases.results")
dat.contr <- read.impute("datasim/OUTPUT_IMPUTE2/imputed.controls.results")
row.names(dat.cases) = sapply(1:dim(dat.cases)[1], function(i) paste("Case",i,sep=""))
row.names(dat.contr) = sapply(1:dim(dat.contr)[1], function(i) paste("Control",i,sep=""))
genotypes = rbind(dat.contr, dat.cases)
phenotypes = c(rep(0,dim(dat.contr)[1]), rep(1,dim(dat.cases)[1]))
#dat <- genotypes
snpsum.col <- col.summary(genotypes)

# thresholds for filtering
call <- 0.95
minor <- 0.01

use <- with(snpsum.col, (!is.na(MAF) & MAF > minor) & Call.rate >= call)
use[is.na(use)] <- FALSE

cat(ncol(genotypes) - sum(use), "SNP's removed due to low MAF or call rate.\n")

genotypes <- genotypes[, use]
snpsum.col <- snpsum.col[use,]

hardy = 10^-6      # HWE cut-off

snpsum.colCont = col.summary(genotypes)
HWEuse = with(snpsum.colCont, !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hardy/2) ) ) )
rm(snpsum.colCont)
HWEuse[is.na(HWEuse)] = FALSE          # Remove NA's as well
cat(ncol(genotypes)-sum(HWEuse),"SNPs will be removed due to high HWE.\n")

genotypes <- genotypes[, HWEuse]
snpsum.col <- snpsum.col[HWEuse,]

X = as(genotypes, "numeric")

Sigma = cov(X)
Sigma.distance = as.dist(1 - abs(cov2cor(Sigma)))
fit = hclust(Sigma.distance, method="single")
corr_max = 0.75
clusters = cutree(fit, h=1-corr_max)

set.seed(123)
ind.screen = rep(F, length(phenotypes))
ind.screen[sample.int(length(phenotypes), size=length(phenotypes)*0.2)] = T
pvals.screen = p.value(single.snp.tests(phenotypes[ind.screen], snp.data = genotypes[ind.screen,]), df=1)

ind.repr = sapply(1:max(clusters), function(c) {
  cluster_elements = clusters==c
  top_within = which.min(pvals.screen[cluster_elements])
  if( length(top_within)==0 ) top_within = 1
  which(cluster_elements)[top_within]
})
X = X[,ind.repr]
pvals.screen = pvals.screen[ind.repr]

library(SNPknock)

Xout_path <- "./datasim/OUTPUT_fastPHASE/X_input"
Xinp_file = SNPknock.fp.writeX(X)
fp_path  = "./datasim/fastPHASE"
fp_output_path = "./datasim/OUTPUT_fastPHASE"
fp_outPath = SNPknock.fp.runFastPhase(fp_path, Xinp_file, out_path = fp_output_path, K=8, numit=10)

r_file = paste(fp_outPath, "_rhat.txt", sep="")
theta_file = paste(fp_outPath, "_thetahat.txt", sep="")
alpha_file = paste(fp_outPath, "_alphahat.txt", sep="")
char_file = paste(fp_outPath, "_origchars", sep="")
hmm = SNPknock.fp.loadFit(r_file, theta_file, alpha_file, X[1,])

Xk = SNPknock.knockoffHMM(X, hmm$pInit, hmm$Q, hmm$pEmit)

plot(colMeans(X),colMeans(Xk),col=rgb(0,0,0,alpha=0.1), pch=16,cex=1); abline(a=0, b=1, col='red', lty=2)

library(knockoff)

Xk[ind.screen,] = X[ind.screen,]
# this yields NAs because some columns (e.g. Xk[, 181] are all 0)


# W = stat.glmnet_coefdiff(X, Xk, phenotypes, family="binomial")
W = stat.stability_selection(X, Xk, phenotypes)
# plot(W, pch=16, cex=1)
t = knockoff.threshold(W, fdr=0.1, offset=0)
discoveries = which(W >= t)
names(discoveries) = colnames(genotypes)[discoveries]
print(discoveries)
colors = rep("gray",length(W))
colors[discoveries] = "blue"
plot(W, col=colors, pch=16, cex=1); abline(h=t, lty=2)

