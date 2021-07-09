#!/bin/bash

#Name the job
#PBS -N StitchNetworkRes

# Options for output and error files
#PBS -j oe

# Flush output file
#PBS -k o

# Use the submission environment
#PBS -V

#Specify walltime, cores, memory
#PBS -l walltime=01:00:00
#PBS -l ncpus=1
#PBS -l mem=4GB

cd  $PBS_O_WORKDIR

echo PBS: node file is $PBS_NODEFILE

# matlab -nodisplay -nosplash -nodesktop -r "combine_output_files('work_pattern_results/infection_output_CT_engagement_withCT',50,20,100000,365,11,'work_pattern_results_combined/infection_output_CT_engagement_withCT.mat');exit"
# matlab -nodisplay -nosplash -nodesktop -r "combine_output_files('work_pattern_results/infection_output_amount_backwards_CT_withCT',50,20,100000,365,11,'work_pattern_results_combined/infection_output_amount_backwards_CT_withCT.mat');exit"
# matlab -nodisplay -nosplash -nodesktop -r "combine_output_files('work_pattern_results/infection_output_run_one_run',25,20,100000,365,1,'work_pattern_results_combined/infection_output_noCT.mat');exit"
# matlab -nodisplay -nosplash -nodesktop -r "combine_output_files('work_pattern_results/infection_output_workplace_CT_threshold_withCT',25,20,100000,365,11,'work_pattern_results_combined/infection_output_workplace_CT_threshold_withCT_var.mat');exit"
matlab -nodisplay -nosplash -nodesktop -r "combine_output_files('work_pattern_results/infection_output_CT_engagement_withCT',25,20,100000,365,11,'work_pattern_results_combined/infection_output_CT_engagement_withCT_var.mat');exit"

echo "Finished"
