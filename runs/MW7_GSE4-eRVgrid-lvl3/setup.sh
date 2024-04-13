#!/bin/bash

Rs=$1
Vv=$2
E=$3

DIR=lvl3-Rs${Rs}-Vv${Vv}-e${E}

if [ -d ${DIR} ]; then
    echo "Directory ${DIR} already exists, exiting..."
    exit 1
fi

mkdir ${DIR}
cd ${DIR}
    echo "doing Rstart=${Rs}, Vvir=${Vv}, e=${E}"

    rsync -ap --exclude arepo ../lvl3-template/* ./
    git clone ../lvl3-template/arepo > git.log 2>&1
    
    cd arepo
    ln -s ../Config.sh ./
    cd ../
    
    sed -i "s/Rs_Vv_e/Rs${Rs}_Vv${Vv}_e${E}/g" job_lvl3.sh job_lvl3-rst.sh
    sed -i "s/Rs-Vv-e/Rs${Rs}-Vv${Vv}-e${E}/g" param_lvl3.txt

    JOBID=`sbatch --parsable job_lvl3.sh`
    sbatch --dependency=afterok:${JOBID} job_lvl3-rst.sh

    echo "done"
cd ../

