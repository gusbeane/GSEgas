module purge
module load gcc
module load openmpi
module load gsl
##module load hdf5
module load python

export LD_LIBRARY_PATH='/n/holylfs05/LABS/hernquist_lab/Users/abeane/hdf5-1.14.1-2/lib':${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH='/n/holylfs05/LABS/hernquist_lab/Users/abeane/code/ffmpeg-6.0/lib':${LD_LIBRARY_PATH}
export PATH='/n/holylfs05/LABS/hernquist_lab/Users/abeane/code/ffmpeg-6.0/bin':${PATH}

export PATH='/n/holylfs05/LABS/hernquist_lab/Users/abeane/code/nasm-2.16.01/bin':$PATH
export LD_LIBRARY_PATH='/n/holylfs05/LABS/hernquist_lab/Users/abeane/code/nasm-2.16.01/lib':$LD_LIBRARY_PATH

export PATH='/n/holylfs05/LABS/hernquist_lab/Users/abeane/code/x264-master/bin':$PATH
export LD_LIBRARY_PATH='/n/holylfs05/LABS/hernquist_lab/Users/abeane/code/x264-master/lib':$LD_LIBRARY_PATH


