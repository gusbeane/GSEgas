#!/bin/sh
##SBATCH -p serial_requeue
#SBATCH -p sapphire,hernquist,conroy,itc_cluster
#SBATCH -J post067

#SBATCH -n 16
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=3800
#SBATCH -t 0-00:20

#SBATCH -o post/slurmlog/OUTPUT.%j_%a.out
#SBATCH -e post/slurmlog/ERROR.%j_%a.err
#SBATCH --open-mode=append

#SBATCH --mail-user=slurmsara@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#SBATCH --array=0-2000:20
##SBATCH --array=0,20

# exit if there are any errors
set -e
set -o pipefail

source ./load-modules.sh

ulimit -c unlimited

START_ID=$((SLURM_ARRAY_TASK_ID))
END_ID=$((SLURM_ARRAY_TASK_ID + 19))

MAX_SNAPNUM=$(ls output/snapshot_[0-9]*.hdf5 | grep -oP 'snapshot_\K[0-9]+' | sort -n | tail -1)
MAX_SNAPNUM=$((10#$MAX_SNAPNUM))

for (( SNAPNUM=START_ID; SNAPNUM<=END_ID; SNAPNUM++ )); do
    if [[ $SNAPNUM -gt $MAX_SNAPNUM ]]; then
        echo "Reached the maximum snap number ($MAX_SNAPNUM), exiting..."
        break
    fi
    ./run_post.sh param_lvl5.txt $SNAPNUM 
done

