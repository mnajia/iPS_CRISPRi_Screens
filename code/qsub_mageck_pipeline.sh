#! /bin/bash
#Mohamad Najia
#Run the MAGeCK pipeline to determine screen hits

current_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )";
repository_dir="${current_dir%/*}"

screen_name=iPS_SSEA5_D9_CRISPRi
mageck_dir=${repository_dir}/SSEA5_differentiation_screen/mageck_output


qsub -N $screen_name -o $mageck_dir/${screen_name}_mageck_pipeline.out.log -e $mageck_dir/${screen_name}_mageck_pipeline.err.log ./mageck_pipeline.sh $screen_name $repository_dir $mageck_dir
