#!/bin/bash

# Milky Way
makenewdiskhalogas/MakeDiskGalaxy MW_lvl5.txt > MND_MW.out 2> MND_MW.err
python3 pass_ics.py MW_ICs.dat MW.hdf5 0

# GSE
makenewdiskhalogas/MakeDiskGalaxy GSE_lvl5.txt > MND_GSE.out 2> MND_GSE.err
python3 pass_ics.py GSE_ICs.dat GSE.hdf5 2

