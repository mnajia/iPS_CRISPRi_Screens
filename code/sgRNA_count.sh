#!/bin/sh
#$ -l h_vmem=15g
#$ -cwd
#$ -q broad
#$ -l h_rt=2:00:00


source /broad/software/scripts/useuse 
reuse .r-3.6.0-bioconductor
export R_LIBS_USER="/broad/blainey_lab/Mo/diffTF_env/R/3.6/libs/"

#parse command line inputs
fq_dir=$1
sgrna_lib_file=$2
out_file=$3

Rscript ./sgRNA_count_v4.1.R $fq_dir $sgrna_lib_file $out_file
