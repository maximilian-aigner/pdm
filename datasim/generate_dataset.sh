#!/bin/bash

n_active_genes=10
n_snps=10000
n_controls=250
n_cases=100

fname_m=genetic_map_chr22_combined_b36.txt
fname_l=1000G/phaseI_b37_no_singletons/ # .legend file
fname_h= # .hap/.haplotype file
output_file=$(basename "$fname_m" .txt)

Rcall=`Rscript --vanilla ../R/determine_simulation_parameters.R $fname_l`

echo $output_file
echo $Rcall
