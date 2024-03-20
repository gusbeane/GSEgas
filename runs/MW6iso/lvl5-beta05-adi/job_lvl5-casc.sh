#!/bin/sh
#SBATCH -p itc_cluster,hernquist,conroy
##SBATCH -p hernquist_ice
#SBATCH -J MW6beta05adi_l5
#SBATCH -n 48
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH -o output/OUTPUT.%j.out
#SBATCH -e output/ERROR.%j.err
#SBATCH --exclusive
##SBATCH --contiguous
#SBATCH --mail-user=slurmsara@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem-per-cpu=3800
#SBATCH -t 7-00:00           # Runtime in D-HH:MM

# exit if there are any errors
set -e
set -o pipefail

source ./load-modules.sh

if [ -f output/hosts ]; then
  rm output/hosts
fi

hostlist=$(scontrol show hostname $SLURM_JOB_NODELIST)

for f in $hostlist
  do
    echo $f':'$SLURM_NTASKS_PER_NODE >> output/hosts
  done

ulimit -c unlimited

cd arepo
make clean > ../output/MAKE.out 2> ../output/MAKE.err
make -j >> ../output/MAKE.out 2>> ../output/MAKE.err
cd ../

mpiexec --mca mpi_leave_pinned 0 --mca oob_tcp_listen_mode listen_thread -np $SLURM_NTASKS  arepo/Arepo param_lvl5.txt 

