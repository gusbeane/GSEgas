#!/bin/bash

Rs=$1
Vv=$2
E=$3

mkdir lvl5-Rs${Rs}-Vv${Vv}-e${E}
cd lvl5-Rs${Rs}-Vv${Vv}-e${E}

    echo "doing Rstart=${Rs}, Vvir=${Vv}, e=${E}"

    rsync -ap --exclude arepo ../lvl5-template/* ./
    git clone ../lvl5-template/arepo > git.log 2>&1
    
    cd arepo
    ln -s ../Config.sh ./
    cd ../
    
    sed -i "s/Rs_Vv_e/Rs${Rs}_Vv${Vv}_e${E}/g" job_lvl5.sh
    sed -i "s/Rs-Vv-e/Rs${Rs}-Vv${Vv}-e${E}/g" param_lvl5.txt

    sbatch job_lvl5.sh

    echo "done"

cd ../

