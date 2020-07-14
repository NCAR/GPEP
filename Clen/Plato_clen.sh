#!/bin/bash

#! Define configuration flags
# See: https://www.acrc.bris.ac.uk/protected/bc4-docs/scheduler/index.html

#SBATCH --job-name=Clen
#SBATCH --time=2-0:0:0
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G

#! add the MATLAB module (as per BCp4)
module load matlab/R2017b
#! Run the job
matlab -nodisplay -r Clen