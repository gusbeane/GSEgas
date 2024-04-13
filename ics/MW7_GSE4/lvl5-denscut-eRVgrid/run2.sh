#!/bin/bash

for Rstart in 116 129 142
do
    for Vvir in 116 129 142
    do
        for e in 04 05 06
        do
            edot="${e:0:1}.${e:1}"
            python3 combine_ics.py MW.hdf5 GSE.hdf5 MW_GSE-Rs${Rstart}-Vv${Vvir}-e${e}.hdf5 ${Rstart} ${Vvir} ${edot}
        done
    done
done

