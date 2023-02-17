#PBS -N buildcalib
#PBS -q share
#PBS -l walltime=1:00:00
#PBS -A P08010000
#PBS -l select=1:ncpus=5

module load conda
module load cdo
module load parallel
conda activate npl-2022b-tgq

python main.py
