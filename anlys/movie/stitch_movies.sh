
#!/bin/bash

ffmpeg -i movies/${1}_star_xz.mp4 -i movies/${1}_gas_xz.mp4 -i movies/${1}_star_xy.mp4 -i movies/${1}_gas_xy.mp4 -filter_complex "[0:v][1:v]hstack[top];[2:v][3:v]hstack[bottom];[top][bottom]vstack[out]" -map "[out]" movies_stitch/${1}.mp4

ffmpeg -i movies/GSE2iso_fg0.7-lvl4-rng_-15_15_-15_15_BoxCenter_nres1024_gas_xz.mp4 -i movies/GSE4iso_fg0.7-lvl4-rng_-15_15_-15_15_BoxCenter_nres1024_gas_xz.mp4 -i movies/GSE5iso_fg0.7-lvl4-rng_-15_15_-15_15_BoxCenter_nres1024_gas_xz.mp4  -filter_complex hstack=inputs=3 movies_stitch/GSE2_GSE4_GSE5.mp4

ffmpeg -i movies/GSE2iso_fg0.7-lvl4-rng_-15_15_-15_15_BoxCenter_nres1024_gas_xz.mp4 -i movies/GSE2iso_fg0.7-lvl4-fSN-rng_-15_15_-15_15_BoxCenter_nres1024_gas_xz.mp4 -i movies/GSE4iso_fg0.7-lvl4-rng_-15_15_-15_15_BoxCenter_nres1024_gas_xz.mp4 -i movies/GSE5iso_fg0.7-lvl4-rng_-15_15_-15_15_BoxCenter_nres1024_gas_xz.mp4  -filter_complex hstack=inputs=4 movies_stitch/GSE2_fSN_GSE4_GSE5.mp4

