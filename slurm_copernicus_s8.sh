#!/bin/bash

#! Define configuration flags
# See: https://www.acrc.bris.ac.uk/protected/bc4-docs/scheduler/index.html

#SBATCH --job-name=oiloo
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --account=hpc_c_giws_clark

#! add the python module
#module load python/3.7.4

# run the application
srun python -u s8_oimerge_200424test.py prcp ${SLURM_ARRAY_TASK_ID}


