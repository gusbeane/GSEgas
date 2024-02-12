#!/bin/sh
##SBATCH -p hernquist,itc_cluster,conroy
#SBATCH -p sapphire
#SBATCH -J gr187_rst
#SBATCH -n 672
#SBATCH -N 6
#SBATCH --ntasks-per-node=112
#SBATCH -o output/OUTPUT.%j.out
#SBATCH -e output/ERROR.%j.err
#SBATCH --mail-user=angus.beane+slurm@cfa.harvard.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem-per-cpu=8600
#SBATCH -t 3-00:00           # Runtime in D-HH:MM

source ./load-modules.sh

ulimit -c unlimited

export SYSTYPE="cannon"

if [ -f "running" ]; then
    echo "Error: 'running' file exists. Exiting."
    exit 1 # Exit with an error code
fi

cd arepo
#make clean > ../output/MAKE.out 2> ../output/MAKE.err
#make -j >> ../output/MAKE.out 2>> ../output/MAKE.err
make_output=$(make -n)
if echo "$make_output" | grep -q "Nothing to be done for 'all'" ; then 
    echo "Code needs no recompiling, continuing"
else
    echo "Code needs to be recompiled, going to exit..."
    echo $make_output
    exit 1
fi

cd ../

mpiexec --mca mpi_leave_pinned 0 --mca oob_tcp_listen_mode listen_thread -np $SLURM_NTASKS  arepo/Arepo param.txt 1 

