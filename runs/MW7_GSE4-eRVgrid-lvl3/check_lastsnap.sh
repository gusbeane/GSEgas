#!/bin/bash

for Rstart in 116 129 142
do
    for Vvir in 116 129 142
    do
        for e in 04 05 06
        do
            ls -lthr lvl3-Rs${Rstart}-Vv${Vvir}-e${e}/output/snap* | tail -1
        done
    done
done

