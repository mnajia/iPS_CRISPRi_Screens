#! /bin/bash
#$ -l h_vmem=20g
#$ -cwd
#$ -q broad
#$ -l h_rt=1:00:00

source /broad/software/scripts/useuse
use Anaconda3

source activate /broad/blainey_lab/Mo/mageck_env


#get script input variables
screen_name=$1
repository_dir=$2
mageck_dir=$3

#declare variables
sgrna_lib_file=${repository_dir}/epifactors_sgRNA_library/EpiFactors_CRISPRi_lib_v3_mageck.csv
control_sgrnas=${repository_dir}/epifactors_sgRNA_library/non_targeting_control_sgRNAs_mageck.txt
fastq_dir=/broad/blainey_lab/Mo/projects/epigenetic_screens/oPS_EpiKDLv3_1_library/iPS_CRISPRi_lib/20230912_1157-2_R1.C03_oPS_EpiKDL_v3_1_SSEA5_Differentiation_Screen_miniseq/fastq


#generate count table
mkdir -p $mageck_dir/count
cd $mageck_dir/count

mageck count \
	--output-prefix $screen_name \
	--norm-method control \
	--list-seq $sgrna_lib_file \
	--fastq $fastq_dir/D9_rep1_SSEA5_hi_S2_R1_001.fastq $fastq_dir/D9_rep1_SSEA5_low_S1_R1_001.fastq $fastq_dir/D9_rep2_SSEA5_hi_S4_R1_001.fastq $fastq_dir/D9_rep2_SSEA5_low_S3_R1_001.fastq $fastq_dir/D9_rep3_SSEA5_hi_S6_R1_001.fastq $fastq_dir/D9_rep3_SSEA5_low_S5_R1_001.fastq \
	--sample-label D9_rep1_SSEA5_hi,D9_rep1_SSEA5_low,D9_rep2_SSEA5_hi,D9_rep2_SSEA5_low,D9_rep3_SSEA5_hi,D9_rep3_SSEA5_low \
	--count-pair False \
	--control-sgrna $control_sgrnas \
	--trim-5 AUTO 


#identify hits with MAGeCK RRA
count_file=${mageck_dir}/count/${screen_name}.count.txt
mkdir -p $mageck_dir/test
cd $mageck_dir/test

mageck test \
	-k $count_file \
	-t D9_rep1_SSEA5_low,D9_rep2_SSEA5_low,D9_rep3_SSEA5_low \
	-c D9_rep1_SSEA5_hi,D9_rep2_SSEA5_hi,D9_rep3_SSEA5_hi \
	-n $screen_name \
	--norm-method none \
	--paired \
	--control-sgrna $control_sgrnas


conda deactivate

