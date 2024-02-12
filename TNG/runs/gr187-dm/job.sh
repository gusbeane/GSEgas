#!/bin/sh
#SBATCH -p hernquist,itc_cluster,conroy
#SBATCH -J gr187-dm
#SBATCH -n 192
#SBATCH -N 4
#SBATCH --ntasks-per-node=48
#SBATCH -o output/OUTPUT.%j.out
#SBATCH -e output/ERROR.%j.err
#SBATCH --mail-user=angus.beane+slurm@cfa.harvard.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem-per-cpu=3800
#SBATCH -t 7-00:00           # Runtime in D-HH:MM

source ./load-modules.sh

ulimit -c unlimited

export SYSTYPE="cannon"

cd arepo
make clean > ../output/MAKE.out 2> ../output/MAKE.err
make -j >> ../output/MAKE.out 2>> ../output/MAKE.err
cd ../

mpiexec --mca mpi_leave_pinned 0 --mca oob_tcp_listen_mode listen_thread -np $SLURM_NTASKS  arepo/Arepo param.txt 

