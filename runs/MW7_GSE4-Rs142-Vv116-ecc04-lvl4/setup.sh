#!/bin/bash

SEED=$1

mkdir lvl4-sd${SEED}
cd lvl4-sd${SEED}

    echo "doing SEED=${SEED}"

    rsync -ap --exclude arepo ../lvl4-template/* ./
    git clone ../lvl4-template/arepo > git.log 2>&1
    
    echo "RANDOM_REALIZATION=${SEED}" >> Config.sh

    cd arepo
    ln -s ../Config.sh ./
    cd ../
    
    sed -i "s/Rs_Vv_e/Rs${Rs}_Vv${Vv}_e${E}/g" job_lvl4.sh

    sbatch job_lvl4.sh

    echo "done"

cd ../

