#!/bin/sh
#SBATCH -p sapphire
#SBATCH -J MW7_GSE4_l3
#SBATCH -n 128 
#SBATCH -N 2
#SBATCH --ntasks-per-node=64
#SBATCH -o output/OUTPUT.%j.out
#SBATCH -e output/ERROR.%j.err
##SBATCH --exclusive
#SBATCH --mail-user=slurmsara@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem-per-cpu=6000
#SBATCH -t 3-00:00           # Runtime in D-HH:MM

source ./load-modules.sh

ulimit -c unlimited

cd arepo
make clean > ../output/MAKE.out 2> ../output/MAKE.err
make -j >> ../output/MAKE.out 2> ../output/MAKE.err
cd ../

mpiexec --mca mpi_leave_pinned 0 --mca oob_tcp_listen_mode listen_thread -np $SLURM_NTASKS  arepo/Arepo param_lvl3.txt 1 
