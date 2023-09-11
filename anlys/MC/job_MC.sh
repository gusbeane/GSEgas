#!/bin/sh
###SBATCH -p itc_cluster,shared,conroy,hernquist,hernquist_ice
#SBATCH -p itc_cluster,shared,conroy,hernquist
#SBATCH -J MC 
#SBATCH -n 48
#SBATCH -N 1
#SBATCH -o logs/OUTPUT.%j.out
#SBATCH -e logs/ERROR.%j.err
#SBATCH --mail-user=angus.beane@cfa.harvard.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem=145G
#SBATCH -t 0-06:00           # Runtime in D-HH:MM

source ../../load-modules.sh

ulimit -c unlimited

python3 compute_MC.py ${SLURM_NTASKS} $1

