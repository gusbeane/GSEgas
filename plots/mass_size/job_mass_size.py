#!/bin/sh
##SBATCH -p itc_cluster,shared,conroy,hernquist,hernquist_ice
#SBATCH -p hernquist,hernquist_ice,itc_cluster,conroy
#SBATCH -J mass_isze 
#SBATCH -n 48
#SBATCH -N 1
#SBATCH -o logs/OUTPUT.%j.out
#SBATCH -e logs/ERROR.%j.err
#SBATCH --mail-user=slurmsara@cfa.harvard.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem=120G
#SBATCH -t 0-06:00           # Runtime in D-HH:MM

source ../../load-modules.sh

ulimit -c unlimited

python3 compute_mass_size.py ${SLURM_NTASKS}

