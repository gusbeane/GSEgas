#!/bin/bash

for sd in {0..15}
do
    ls -lthr lvl4-sd${sd}/output/snap* | tail -1
done

