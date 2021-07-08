#!/bin/bash
#SBATCH --job-name=COVIDmodel
#SBATCH --time=01:00:00
#SBATCH --mem=4g
#SBATCH --ntasks=1
#SBATCH --export=ALL
#SBATCH --partition=ntd


# ${SLURM_SUBMIT_DIR} points to the path where this script was submitted from
cd ${SLURM_SUBMIT_DIR}

module load matlab

matlab -nodisplay -nosplash -nodesktop -r "combine_output_files_all;exit"

echo "Finished"
