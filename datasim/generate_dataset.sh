#!/bin/bash

n_active_genes=10
n_snps=99999
n_controls=100
n_cases=100

# set up various paths

# ----
# EDIT THESE:
files_dir=hapmap3_r2_b36
fname_m=$files_dir/genetic_map_chr22_combined_b36.txt
imputed_name=hapmap3_r2_b36_chr22
# ----

# These should not need editing, but sometimes the extension is different
fname_l=$files_dir/$imputed_name.legend # .legend file
fname_h=$files_dir/$imputed_name.haps # .hap/.haplotype file

output_dir=working_dataset
output_file=$(basename "$fname_m" .txt)

Rcall_args=`Rscript --vanilla ../R/determine_simulation_parameters.R $fname_l $n_snps`
echo $Rcall_args

# build HAPGEN2 command string
hapgen2_call_str="./hapgen2 -m $fname_m -l $fname_l -h $fname_h -o $output_dir/hapgen2/$output_file -n $n_controls $n_cases $Rcall_args"
echo $hapgen2_call_str
hapgen2_call=eval $hapgen2_call_str

# prepare to run IMPUTE2 - twice
# try to run first with -dl and hope IMPUTE2 silently ignores it
just_intervals="-$(cut -d'-' -f2 <<<$Rcall_args)"
echo $just_intervals

impute2_call_str="./impute2 -h $fname_h -l $fname_l -m $fname_m -g $output_dir/hapgen2/$output_file.controls.gen -o $output_dir/impute2/$output_file.imputed.controls $just_intervals"
echo $impute2_call_str
impute2_call=eval $impute2_call_str

impute2_call_str2="./impute2 -h $fname_h -l $fname_l -m $fname_m -g $output_dir/hapgen2/$output_file.cases.gen -o $output_dir/impute2/$output_file.imputed.cases $just_intervals"
echo $impute2_call_str2
impute2_call2=eval $impute2_call_str2

