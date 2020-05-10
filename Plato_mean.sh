#!/bin/bash

#! Define configuration flags
# See: https://www.acrc.bris.ac.uk/protected/bc4-docs/scheduler/index.html

#SBATCH --job-name=minmax
#SBATCH --time=0-6:00:00
#SBATCH --mem-per-cpu=10G

#! add the python module
module load python/3.7.4

# run the application
srun python -u minmax_to_meanrange.py




