#!/bin/bash

for Rstart in 116 129 142
do
    for Vvir in 116 129 142
    do
        for e in 04 05 06
        do
            ./setup.sh ${Rstart} ${Vvir} ${e}
        done
    done
done

