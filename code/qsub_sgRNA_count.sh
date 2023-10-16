#! /bin/bash
#Mohamad Najia
#Quantify abundance of sgRNAs from NGS reads

current_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )";
repository_dir="${current_dir%/*}"

sgrna_lib_file=${repository_dir}/epifactors_sgRNA_library/EpiFactors_CRISPRi_lib_v3.tsv
samples_file=${repository_dir}/SSEA5_differentiation_screen/samples.txt
output_dir=${repository_dir}/SSEA5_differentiation_screen/sgRNA_counts


while IFS=$'\t' read  sample  fq_file
do
  echo "Processing sample: " $sample 
  output_file=$output_dir/${sample}_counts.rds
  qsub -N $sample -o $output_dir/${sample}_count.out.log -e $output_dir/${sample}_count.err.log ./sgRNA_count.sh $fq_file $sgrna_lib_file $output_file
done < ${samples_file}

