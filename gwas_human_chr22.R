library(snpStats)
set.seed(43192)
library(plyr)


dat.cases <- read.impute("datasim/working_dataset/hapgen2/generated_output.cases.gen")
# dat.cases <- read.impute("datasim/tutorial_data/sim.out.cases.gen")
dat.contr <- read.impute("datasim/working_dataset/hapgen2/generated_output.controls.gen")
# dat.contr <- read.impute("datasim/tutorial_data/sim.out.controls.gen")
row.names(dat.cases) = sapply(1:dim(dat.cases)[1], function(i) paste("Case",i,sep=""))
row.names(dat.contr) = sapply(1:dim(dat.contr)[1], function(i) paste("Control",i,sep=""))
genotypes = rbind(dat.contr, dat.cases)
phenotypes = c(rep(0,dim(dat.contr)[1]), rep(1,dim(dat.cases)[1]))
#dat <- genotypes
snpsum.col <- col.summary(genotypes)

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

# annotation step
legend.file <- read.csv("datasim/working_dataset/hapgen2/generated_output.legend", header = TRUE, sep = ' ', stringsAsFactors = FALSE)
legend.file <- legend.file[use, ]
legend.file <- legend.file[HWEuse, ]
ann.df <- data.frame(rsid=legend.file$rs, chr=rep("chr22", nrow(legend.file)), pos = legend.file$pos)

source('annotation.R')
annotated <- annotate_snps_with_genes(ann.df)

#baseann <- cbind(ann.df, ann.df$rsid)
#colnames(baseann) <- c("rsid", "chr", "pos", "grp")

#just.groups <- data.frame(rsid=annotated$names,
#                           grp=annotated$GENESYMBOL)
#just.groups <- data.frame(lapply(just.groups, as.character), stringsAsFactors=FALSE)
#just.groups[is.na(just.groups$grp), ]$grp <- just.groups[is.na(just.groups$grp),]$rsid

#avec <- just.groups$grp
#names(avec)<-just.groups$rsid
annotated <- annotated[complete.cases(annotated),]


group.names <- colnames(X)
# not_subbed <- group.names[!(group.names %in% annotated$names)]
group.names <- mapvalues(group.names, from=annotated$names, to=annotated$GENESYMBOL)
# map <- setNames(c(annotated$names, not_subbed), c(annotated$GENESYMBOL, not_subbed))
# group.names[] <- map[unlist(group.names)]
# idx <- annotated$names == group.names
# group.names[idx] <- annotated$GENESYMBOL[idx]
#group.names[annotated$names %in% group.names] <- annotated[annotated$names %in% group.names, ]$GENESYMBOL
  
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

plot(colMeans(X),colMeans(Xk),col = rgb(0,0,0,alpha = 0.1), pch=16,cex=1); abline(a=0, b=1, col='red', lty=2)

library(knockoff)
library(grpreg)

total.groups <- c(group.names, paste(group.names, "_knockoff"))
grp.fit <- cv.grpreg(cbind(X, Xk), phenotypes, total.groups, family = "binomial", penalty="cMCP")
