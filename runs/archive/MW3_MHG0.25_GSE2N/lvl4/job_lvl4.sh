#!/bin/sh
##SBATCH -p itc_cluster,shared,conroy
#SBATCH -p hernquist_ice
#SBATCH -J MW3GSE2Nmerge0l4
#SBATCH -n 256
#SBATCH -N 4
#SBATCH --ntasks-per-node=64
#SBATCH -o output/OUTPUT.%j.out
#SBATCH -e output/ERROR.%j.err
#SBATCH --exclusive
##SBATCH --contiguous
#SBATCH --mail-user=angus.beane@cfa.harvard.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem-per-cpu=7900
#SBATCH -t 7-00:00           # Runtime in D-HH:MM

source ./load-modules.sh

if [ -f hosts ]; then
  rm hosts
fi

hostlist=$(scontrol show hostname $SLURM_JOB_NODELIST)

for f in $hostlist
  do
    echo $f':'$SLURM_NTASKS_PER_NODE >> hosts
  done

ulimit -c unlimited

cd arepo
make clean ; make -j
cd ../

mpiexec --mca mpi_leave_pinned 0 --mca oob_tcp_listen_mode listen_thread -np $SLURM_NTASKS  arepo/Arepo param_lvl4.txt 

