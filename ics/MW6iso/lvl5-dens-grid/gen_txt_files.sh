#!/bin/bash

for beta in 08 07 067 06 05
do
    mkdir -p beta${beta}
    cd beta${beta}

    for RC in 5 9 20 30 40
    do
        for vphi in 01 02 03
        do
            mkdir -p RC${RC}-vphi${vphi}
            cd RC${RC}-vphi${vphi}

            cp ../../MW_lvl5-beta${beta}.txt MW_lvl5-beta${beta}-RC${RC}-vphi${vphi}.txt
            
            vphi_num="${vphi:0:1}.${vphi:1}"
            sed -i "s/VelPhiFactor            0.2/VelPhiFactor            ${vphi_num}/g" MW_lvl5-beta${beta}-RC${RC}-vphi${vphi}.txt
            sed -i "s/RC                      9.0/RC                      ${RC}/g" MW_lvl5-beta${beta}-RC${RC}-vphi${vphi}.txt
            sed -i "s/OutputFile              MW_ICs-beta${beta}.dat/OutputFile              MW_ICs-beta${beta}-RC${RC}-vphi${vphi}.dat/g" MW_lvl5-beta${beta}-RC${RC}-vphi${vphi}.txt

            cd ../
        done
    done

    cd ../
done

