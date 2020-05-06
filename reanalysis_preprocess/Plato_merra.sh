#!/bin/bash

#! Define configuration flags
# See: https://www.acrc.bris.ac.uk/protected/bc4-docs/scheduler/index.html

#SBATCH --job-name=merrapre
#SBATCH --cpus-per-task=6
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=20G

#! add the python module
module load matlab/R2017b

# run the application
matlab -r -nodisplay MERRA2_read_plato