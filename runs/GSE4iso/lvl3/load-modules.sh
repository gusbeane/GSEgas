module purge
module load gcc
module load openmpi
module load gsl
#module load hdf5
module load python

export LD_LIBRARY_PATH='/n/holylfs05/LABS/hernquist_lab/Users/abeane/hdf5-1.14.1-2/lib':${LD_LIBRARY_PATH}

#source /n/holylfs05/LABS/hernquist_lab/Lab/spack/share/spack/setup-env.sh

#spack load gcc
#spack load openmpi
#spack load gsl
#spack load hdf5
#spack load python

