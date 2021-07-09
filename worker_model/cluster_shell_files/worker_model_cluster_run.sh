#!/bin/bash
#SBATCH --job-name=uni_model
#SBATCH --time=024:00:00
#SBATCH --mem=2g
#SBATCH --export=ALL
#SBATCH --partition=ntd

# Submit a job array
#SBATCH --array=1-50
#SBATCH --ntasks=1

# ${SLURM_SUBMIT_DIR} points to the path where this script was submitted from
cd ${SLURM_SUBMIT_DIR}

module load julia/1.5.1

# args/ARGS list
# args[1] JOB_ID
# args[2] RNGseed: To be used to initialise the random number generator
# args[3] COUNTFINAL: Number of replicates requested
# args[4] ENDTIME: Timesteps for each individual simulation
# args[5] n_nodes: Overall size of network
# args[6] n_sectors: Defining total number of work sectors
# args[7] runsets: Scenarios to run: [configuration, intervention] (see load_configurations() and load_interventions() for options)
julia worker_model.jl ${SLURM_ARRAY_TASK_ID} 1 100 10 365 10000 41 '["default" "none"]'

module load matlab

matlab -nodisplay -nosplash -nodesktop -r "resave_MAT_file('../../Results/worker_model/',${SLURM_ARRAY_TASK_ID});exit"

echo "Finished"
