#!/bin/sh
#SBATCH -p itc_cluster,shared,conroy,hernquist
#SBATCH -J movie 
#SBATCH -n 48
#SBATCH -N 1
#SBATCH -o logs/OUTPUT_frames.%j.out
#SBATCH -e logs/ERROR_frames.%j.err
#SBATCH --mail-user=angus.beane@cfa.harvard.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem=120G
#SBATCH -t 0-06:00           # Runtime in D-HH:MM

source ../../load-modules.sh

ulimit -c unlimited

python3 compute_movie.py ${SLURM_NTASKS} $1

