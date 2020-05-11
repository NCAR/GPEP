#!/bin/bash

#! Define configuration flags
# See: https://www.acrc.bris.ac.uk/protected/bc4-docs/scheduler/index.html

#SBATCH --job-name=PyGMET
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=20G

#! add the python module
module load python/3.7.4

# run the application
srun python main_CAI.py 19790101 19790131




