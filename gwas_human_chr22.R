library(snpStats)
set.seed(43192)
library(plyr)


dat.cases <- read.impute("datasim/working_dataset/hapgen2/generated_output.cases.gen")
dat.contr <- read.impute("datasim/working_dataset/hapgen2/generated_output.controls.gen")
row.names(dat.cases) = sapply(1:dim(dat.cases)[1], function(i) paste("Case",i,sep=""))
row.names(dat.contr) = sapply(1:dim(dat.contr)[1], function(i) paste("Control",i,sep=""))
genotypes = rbind(dat.contr, dat.cases)
phenotypes = c(rep(0,dim(dat.contr)[1]), rep(1,dim(dat.cases)[1]))
snpsum.col <- col.summary(genotypes)

call <- 0.95
minor <- 0.01

use <- with(snpsum.col, (!is.na(MAF) & MAF > minor) & Call.rate >= call)
use[is.na(use)] <- FALSE

cat(ncol(genotypes) - sum(use), "SNP's removed due to low MAF or call rate.\n")

genotypes <- genotypes[, use]
snpsum.col <- snpsum.col[use,]

hardy = 10^-6

snpsum.colCont = col.summary(genotypes)
HWEuse = with(snpsum.colCont, !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hardy/2) ) ) )
rm(snpsum.colCont)
HWEuse[is.na(HWEuse)] = FALSE          # Remove NA's as well
cat(ncol(genotypes)-sum(HWEuse),"SNPs will be removed due to high HWE.\n")

genotypes <- genotypes[, HWEuse]
snpsum.col <- snpsum.col[HWEuse,]

X = as(genotypes, "numeric")

# Use clusters as groups
Sigma = cov(X)
Sigma.distance = as.dist(1 - abs(cov2cor(Sigma)))
fit = hclust(Sigma.distance, method="single")
corr_max = 0.75
clusters = cutree(fit, h=1-corr_max)

library(pheatmap)
pdf(file = "heatmap.pdf", width = 5, height = 10)
pheatmap(cov2cor(Sigma), cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

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

plot(colMeans(X),colMeans(Xk),col = rgb(0,0,0,alpha = 0.1), pch=16,cex=1);
abline(a=0, b=1, col='red', lty=2)

library(knockoff)
library(grpreg)
library(glmnet)

total.groups <- c(clusters, paste0(clusters, "_knockoff"))
names(total.groups) <- c(names(clusters), paste0(names(clusters), "_knockoff"))
#grp.fit <- grpreg(cbind(X, Xk), phenotypes, total.groups, family = "binomial",
#                  penalty="grSCAD", nlambda = 100)
grp.fit <- glmnet(cbind(X, Xk), phenotypes, family = "binomial", nlambda = 100)
grp.lambdas <- grp.fit$lambda
min_coefs <- 20
n_nzcoefs <- sapply(grp.lambdas, function(l) sum(coef(grp.fit, s = l)!=0))
chosen.lambda <- max(grp.lambdas[n_nzcoefs>=min_coefs])
# chosen.lambda <- grp.fit$lambda.min
Z = coef(grp.fit, s = chosen.lambda)
p = dim(X)[2]
colSD <- apply(cbind(X, Xk), 2, sd)
orig = 2:(p+1)
# Z <- Z[2:length(Z)]
# Z <- Z * c(1,colSD);
W = abs(Z[orig]) - abs(Z[p+orig])

t = knockoff.threshold(W, fdr = 0.1, offset = 0)
discoveries = which(W >= t)
names(discoveries) = colnames(genotypes)[discoveries]
print(discoveries)

real <- read.table("datasim/working_dataset/active_genes.txt", header = TRUE, stringsAsFactors = FALSE)
signals <- real$rs
signals.id <- match(signals, colnames(X))
names(signals.id) <- signals

# Plot W-statistic
colors = rep("gray",length(W))
colors[discoveries] = "red"
colors[signals.id] = "green"
plot(W, col=colors, pch=16, cex=1); abline(h=t, lty=2)

