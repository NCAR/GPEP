#!/bin/bash

#! Define configuration flags
# See: https://www.acrc.bris.ac.uk/protected/bc4-docs/scheduler/index.html

#SBATCH --job-name=PyGMET
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-2:00:00
#SBATCH --mem-per-cpu=10G

#! add the python module
module load python/3.7.4

# run the application
srun python -u s2_find_nearstn.py




