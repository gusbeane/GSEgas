#!/bin/bash

export AREPO_commit=fd85d560d08d5dc76d57733e3f54bfe1f4ea56d1
export REPLACE=fidCOM
export ICS_DIR=../../../ics/fgGSE0.5_fgMW0.5-COM

source load-modules.sh

#for i in 5 4 3 2 1
for i in 4 3 
do
    mkdir lvl${i}
    cd lvl${i}

    cp ../load-modules.sh .

    ln -s ${ICS_DIR}/lvl${i} ICs

    # clone arepo
    cp ../param/param_lvl${i}.txt ./
    git clone git@bitbucket.org:volkerspringel/arepo.git

    cd arepo
    git checkout ${AREPO_commit}
    cp ../../Config-SMUGGLE.sh Config.sh
    cp ../../Makefile.systype .
    # make -j # commented out bc will be recompiled at runtime anyways
    cd ../
    
    cp ../job_lvlN.sh job_lvl${i}.sh
    sed -i "s/REPLACE-/${REPLACE}/g" job_lvl${i}.sh
    sed -i "s/REPli/l${i}/g" job_lvl${i}.sh
    sed -i "s/lvlN/lvl${i}/g" job_lvl${i}.sh

    #mkdir ${SCRATCH}/hernquist_lab/abeane/output-${REPLACE}-l${i}
    #ln -s ${SCRATCH}/hernquist_lab/abeane/output-${REPLACE}-l${i} output
    mkdir output

    cd ../

    echo "done with lvl${i}"
done

