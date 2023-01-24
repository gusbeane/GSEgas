#!/bin/bash

if (( $# != 1 ))
then
  echo "Usage: ./gen_ics.sh i"
  echo "Where i is the desired lvl"
  exit 1
fi

export LVL=lvl$1

source load-modules.sh

MakeNewDisk/MakeDiskGalaxy MW_${LVL}.txt 2>&1 | tee MND_output.txt

MakeNewDisk-no_gas/MakeDiskGalaxy Sgr_${LVL}.txt 2>&1 | tee MND-Sgr_output.txt

mpirun -np 1 arepo_BG/Arepo param_ICs.txt 2>&1 | tee BG_output.txt

python3 create_ics.py 2>&1 | tee create_ics_output.txt

