#!/bin/bash

ICS_DIR=/n/holylfs05/LABS/hernquist_lab/Users/abeane/GSEgas/ics/MW4_MHG0.25_vphi0.2_GSE2_MHG0.5/lvl4-vary
COMMIT=763e6294d
BRANCH=SMUGGLEdev-dt_AB
CONFIG=Config.sh
LVL=4

source load-modules.sh

for icpath in `ls ${ICS_DIR}/ics_*`
do
    icname=$(basename $icpath)
    icname=${icname#ics_}
    icname=${icname%.hdf5}
    echo $icname

    mkdir lvl${LVL}-$icname
    cd lvl${LVL}-$icname

        cp ../load-modules.sh .
   
        ln -s $icpath ics.hdf5
    
        cp ../param/param_lvl${LVL}.txt ./
        cp ../$CONFIG Config.sh

        git clone --depth=32 --branch=${BRANCH} git@bitbucket.org:volkerspringel/arepo.git
    
        cd arepo
            git checkout $COMMIT || exit 1
            ln -s ../Config.sh ./
        cd ../

    for jn in sapph ice casc
    do
        cp ../job-${jn}.sh ./ 
        sed -i "s/REPLACE-/${icname}/g" job-${jn}.sh
        sed -i "s/REPli/l${LVL}/g" job-${jn}.sh
        sed -i "s/lvlN/lvl${LVL}/g" job-${jn}.sh
    done
    
    mkdir output
    echo "done with lvl${LVL}-${icname}"

    cd ../
done


    #mkdir output

    #cd ../

    #echo "done with lvl${i}"



