#!/bin/bash

Rvir_fid=129
Vvir_fid=129
e_fid=0.5
pro_fid=1
angle_fid=-165

# vary Vvir
for Vvir in 114 129 144
do
    for e in 0.3 0.5 0.6
    do
        Rvir=${Rvir_fid}
        pro=${pro_fid}
        angle=${angle_fid}
        outname=ics_Rvir${Rvir}_Vvir${Vvir}_e${e}_pro${pro}_angle${angle}.hdf5

        python3 create_ics.py ${Rvir} ${Vvir} ${e} ${pro} ${angle} ${outname}
    done
done

