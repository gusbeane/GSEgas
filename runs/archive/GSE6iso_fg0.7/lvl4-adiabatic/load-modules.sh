module purge
module load gcc/12.2.0-fasrc01 
module load openmpi
module load gsl
##module load hdf5
module load python

export LD_LIBRARY_PATH='/n/holylfs05/LABS/hernquist_lab/Users/abeane/hdf5-1.14.1-2/lib':${LD_LIBRARY_PATH}

