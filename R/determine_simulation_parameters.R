#!/usr/bin/env Rscript
set.seed(43193)
args = commandArgs(trailingOnly = TRUE)
# choose 'causal' gene loci
k = 5
# N = 100000

legend_file <- args[1]
N <- args[2]
loci.df <- read.csv(paste0('~/src/pdm/datasim/', legend_file), header = TRUE, sep = ' ')

bounded <- head(loci.df[order(loci.df$position), ], n = N)
annotations <- read.table('~/src/pdm/datasim/working_dataset/annotations.txt', stringsAsFactors = FALSE, header = TRUE, sep = ' ')
merged <- merge(bounded, annotations, by.x = "rs", by.y = "names")
# sample two genes, and then k among SNPs
genes <- table(merged$GENESYMBOL)
genes <- genes[genes>k]
top2 <- sample(1:length(genes), 2, replace = FALSE)
selected.genes <- genes[top2]
selected.genes.snps <- subset(merged, GENESYMBOL %in% names(selected.genes))
# merged.2 has at least 20 SNPs

causal.genes.snps <- sample(1:nrow(selected.genes.snps), k, replace = FALSE) # bounded, k, replace = FALSE)
final.selection <- selected.genes.snps[causal.genes.snps, ]
causal.genes.loci <- final.selection$position
loc.interval <- c(min(bounded$position), max(bounded$position)) # might be narrow than (lower, upper)
heter.sizes <- round(runif(k, min = 1.25, max = 1.5), 4)
homoz.sizes <- heter.sizes + round(runif(k, min = 1.25, max = 1.5), 4)
risk.alleles <- rep(0, k)
# risk.alleles <- sample(c(0, 1), replace = TRUE, size = k)

ints <- paste0("-int ", paste(loc.interval, sep = " ", collapse = " "))
dls <- paste("-dl", paste(causal.genes.loci, risk.alleles, heter.sizes, homoz.sizes, sep = " ", collapse = " "))
cat(paste(ints, dls, sep = " ", collapse = " "))

write.table(cbind(final.selection, heter.sizes, homoz.sizes, risk.alleles), "~/src/pdm/datasim/working_dataset/active_genes.txt", quote = FALSE, row.names = FALSE)

# ./hapgen2 -m genetic_map_chr22_combined_b36.txt -l CEU.0908.chr22.legend \
# -h CEU.0908.chr22.hap -o OUTPUT_HAPGEN2/chr22.out -n 50 50 \
# -int 14431249 49588215 \
# -dl 16042533 1 1.4386 1.4163 29977662 1 1.4073 1.2737 42519023 1 1.4275 1.346 18716857 1 1.2502 1.3186 34152478 1 1.3688 1.4537 21692266 1 1.305 1.3621 18914875 1 1.345 1.4525 40915268 1 1.4032 1.4531 45898329 1 1.3379 1.4486 27596224 1 1.2778 1.36
# POTENTIALLY: -t [file containing the set wanted SNP locations, one per line]

# Now, use IMPUTE (DO WE EVEN NEED THIS? — YES)
# ../impute2 -h CEU.0908.chr22.hap -l CEU.0908.chr22.legend -m genetic_map_chr22_combined_b36.txt -g OUTPUT_HAPGEN2/chr22.out.controls.gen -int 20000000 25000000 -o OUTPUT_IMPUTE2/imputed.controls.results
