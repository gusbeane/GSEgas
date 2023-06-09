module purge
module load gcc/7.1.0-fasrc01
module load openmpi/3.1.3-fasrc01
#module load gcc/12.1.0-fasrc01
#module load openmpi/4.1.3-fasrc01
module load gsl
module load hdf5
module load python/3.6.3-fasrc02
module load ffmpeg/4.0.2-fasrc01

export HWLOC_INCL='-I/n/home01/abeane/local/hwloc/include'
export HWLOC_LIB='-L/n/home01/abeane/local/hwloc/lib'
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/n/home01/abeane/local/hwloc/lib

