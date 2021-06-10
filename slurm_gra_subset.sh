#!/bin/bash
#SBATCH --job-name=subset
#SBATCH --time=0-5:00:00
#SBATCH --mem=10G
#SBATCH --account=rpp-kshook
module load python/3.7.4
source ~/ENV/bin/activate

ID=${SLURM_ARRAY_TASK_ID}
step=5
y1=$(((ID-1)*step+1979))
y2=$((y1+step-1))

srun python -u s12_subsetting.py $y1 $y2
