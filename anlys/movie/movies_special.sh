#!/bin/bash

# MW2iso_fg0.5-lvl4-rng_-8_8_-8_8_gas_xy_novSFE1.mp4
ffmpeg -i movies/MW2iso_fg0.5-lvl4-rng_-8_8_-8_8_gas_xy.mp4 -i movies/MW2iso_fg0.5-lvl4-nov-rng_-8_8_-8_8_gas_xy.mp4 -i movies/MW2iso_fg0.5-lvl4-SFE1-rng_-8_8_-8_8_gas_xy.mp4 -i movies/MW2iso_fg0.5-lvl4-novSFE1-rng_-8_8_-8_8_gas_xy.mp4 -filter_complex "[0:v][1:v]hstack[top];[2:v][3:v]hstack[bottom];[top][bottom]vstack[out]" -map "[out]" movies_special/MW2iso_fg0.5-lvl4-rng_-8_8_-8_8_gas_xy_novSFE1.mp4
