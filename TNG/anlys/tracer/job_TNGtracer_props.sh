#!/bin/sh
##SBATCH -p itc_cluster,shared,conroy,hernquist,hernquist_ice
#SBATCH -p hernquist,hernquist_ice,itc_cluster,conroy
#SBATCH -J TNGtrace 
#SBATCH -n 48
#SBATCH -N 1
#SBATCH -o logs/OUTPUT.%j.out
#SBATCH -e logs/ERROR.%j.err
#SBATCH --mail-user=angus.beane@cfa.harvard.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem=120G
#SBATCH -t 0-06:00           # Runtime in D-HH:MM
#SBATCH --array=0-99

source ../../load-modules.sh

ulimit -c unlimited

echo ${SLURM_ARRAY_TASK_ID}
python3 compute_TNGtracer.py ${SLURM_NTASKS} $1 ${SLURM_ARRAY_TASK_ID}

