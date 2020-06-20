#!/bin/bash

#! Define configuration flags
# See: https://www.acrc.bris.ac.uk/protected/bc4-docs/scheduler/index.html

#SBATCH --job-name=corrmerge
#SBATCH --time=0-12:00:00
#SBATCH --mem-per-cpu=20G

#! add the python module
module load python/3.7.4

# run the application
srun python -u temprun.py trange BMA Add_Climo




