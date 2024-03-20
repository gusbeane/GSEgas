#!/bin/sh
#SBATCH -p itc_cluster,conroy,hernquist,hernquist_ice
##SBATCH -p hernquist_ice
#SBATCH -J movie 
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --array=0-1000
#SBATCH -o logs/OUTPUT_frames.%j.out
#SBATCH -e logs/ERROR_frames.%j.err
#SBATCH --mail-user=slurmsara@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem=2G
#SBATCH -t 0-00:15           # Runtime in D-HH:MM

source ../../load-modules.sh

ulimit -c unlimited

python3 compute_movie.py $1 $2 ${SLURM_ARRAY_TASK_ID}

