#!/bin/sh
#SBATCH -p itc_cluster,shared,conroy,hernquist,hernquist_ice
##SBATCH -p hernquist_ice
#SBATCH -J movie 
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --array=0-1000
#SBATCH -o logs/OUTPUT_frames.%j.out
#SBATCH -e logs/ERROR_frames.%j.err
#SBATCH --mail-user=angus.beane@cfa.harvard.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem=2G
#SBATCH -t 1-00:00           # Runtime in D-HH:MM

source ../../load-modules.sh

ulimit -c unlimited

python3 compute_movie.py $1 $2 ${SLURM_ARRAY_TASK_ID}

