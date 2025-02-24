module purge
module load gcc/13.2.0-fasrc01
module load openmpi/5.0.2-fasrc02
module load gsl/2.8-fasrc01
##module load hdf5
module load python/3.12.5-fasrc01

export LD_LIBRARY_PATH='/n/holylfs05/LABS/hernquist_lab/Users/abeane/hdf5-1.14.1-2/lib':${LD_LIBRARY_PATH}

