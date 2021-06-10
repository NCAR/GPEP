#!/bin/bash
#SBATCH --job-name=stnreg
#SBATCH --time=1-0:00:00
#SBATCH --mem=30G
#SBATCH --account=rpp-kshook

module load python/3.7.4
source ~/ENV/bin/activate
srun python -u s3_stn_regression_newpredictor.py $((${SLURM_ARRAY_TASK_ID}+1978))
