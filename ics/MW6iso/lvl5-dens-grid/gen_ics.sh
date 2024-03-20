#!/bin/bash

if (( $# != 1 )); then
    echo -e "Need beta value" >&2
    exit 1  # Return with an error status
fi

beta=$1

cd beta${beta}

for RC in 5 9 20 30 40
do
    for vphi in 01 02 03
    do
        cd RC${RC}-vphi${vphi}
            echo "working on beta=${beta}, RC=${RC}, vphi=${vphi}"
            #../../makenewdiskhalogas/MakeDiskGalaxy MW_lvl5-beta${beta}-RC${RC}-vphi${vphi}.txt > MND.out 2> MND.err 
            python3 ../../pass_MW_ics.py 5 MW_ICs-beta${beta}-RC${RC}-vphi${vphi}.dat MW-beta${beta}-RC${RC}-vphi${vphi}.hdf5
            echo "done"
        cd ../
    done
done

cd ../

