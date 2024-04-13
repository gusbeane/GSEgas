#!/bin/sh
#SBATCH -p itc_cluster,hernquist,conroy,hernquist_ice,sapphire
#SBATCH -J MW7_GSE4_Rs142_Vv129_e06_l5
#SBATCH -n 32
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -o output/OUTPUT.%j.out
#SBATCH -e output/ERROR.%j.err
#SBATCH --mail-user=slurmsara+grid@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem-per-cpu=3800
#SBATCH -t 3-00:00           # Runtime in D-HH:MM

source ./load-modules.sh

ulimit -c unlimited

cd arepo
make clean > ../output/MAKE.out 2> ../output/MAKE.err
make -j >> ../output/MAKE.out 2> ../output/MAKE.err
cd ../

mpiexec --mca mpi_leave_pinned 0 --mca oob_tcp_listen_mode listen_thread -np $SLURM_NTASKS  arepo/Arepo param_lvl5.txt 

