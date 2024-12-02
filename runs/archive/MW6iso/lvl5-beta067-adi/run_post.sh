#!/bin/bash

# cd arepo ; make -j CONFIG=Config_post.sh BUILD=build_post EXEC=Arepo_post ; cd ../

# Step 1: Check for at least two arguments
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <PARAM> <SNAPNUM> [NCPUS]"
    exit 1
fi

# Assign the first two parameters
PARAM="$1"
SNAPNUM="$2"

# If a third parameter is provided, use it as NCPUS; otherwise, default to SLURM_NTASKS
if [ -n "$3" ]; then
    NCPUS="$3"
else
    NCPUS="${SLURM_NTASKS}"
fi

# Step 2: Make directories 
mkdir -p post/snapbak
mkdir -p post/log
mkdir -p 'post/done'

# Step 3: Process snapnum and make sure snapshot exists. Exit if this snap has already been processed
printf -v PADDED_SNAPNUM "%03d" "$SNAPNUM"
SNAPNAME="output/snapshot_${PADDED_SNAPNUM}.hdf5"

if ! test -f $SNAPNAME; then
  echo "Could not find snapshot at $SNAPNAME exiting..."
  exit 1
fi

if test -f "post/done/${PADDED_SNAPNUM}"; then
  echo "Already did snapshot ${PADDED_SNAPNUM}, exiting..."
  exit 0
fi

# Step 4: Make a backup and run postprocessing
cp -p "$SNAPNAME" "post/snapbak/"
    
echo "Working on snapshot number ${PADDED_SNAPNUM}..."

if ! mpirun -np "${NCPUS}" arepo/Arepo_post "$PARAM" 18 "$SNAPNUM" > "post/log/${PADDED_SNAPNUM}.log" 2> "post/log/${PADDED_SNAPNUM}.err"; then
    echo "mpirun command failed for snapshot $PADDED_SNAPNUM. Exiting."
    exit 1
fi

# Step 5: Replace the old snapshot
# mv "output/snapshot_potupdated_${PADDED_SNAPNUM}.hdf5" "$SNAPNAME"

# Step 6: Create done file and exit
touch post/done/${PADDED_SNAPNUM}

echo "Processing complete." 

