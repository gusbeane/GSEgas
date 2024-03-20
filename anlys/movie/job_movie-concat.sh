#!/bin/sh
#SBATCH -p itc_cluster,conroy,hernquist,hernquist_ice
##SBATCH -p hernquist_ice
#SBATCH -J movie 
#SBATCH -n 6
#SBATCH -N 1
#SBATCH -o logs/OUTPUT_frames.%j.out
#SBATCH -e logs/ERROR_frames.%j.err
#SBATCH --mail-user=slurmsara@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem=180G
#SBATCH -t 0-01:00           # Runtime in D-HH:MM

source ../../load-modules.sh

ulimit -c unlimited

python3 compute_movie.py $1 $2 -1

