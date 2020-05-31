#!/bin/bash

#! Define configuration flags
# See: https://www.acrc.bris.ac.uk/protected/bc4-docs/scheduler/index.html

#SBATCH --job-name=readown
#SBATCH --time=0-2:00:00
#SBATCH --mem-per-cpu=10G

#! add the python module
module load python/3.7.4

# run the application
srun python -u reanalysis_downscale2.py $((${SLURM_ARRAY_TASK_ID}+1979)) $((${SLURM_ARRAY_TASK_ID}+1979)) GWR




