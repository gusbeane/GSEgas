
#!/bin/bash

ffmpeg -i movies/${1}_star_xz.mp4 -i movies/${1}_gas_xz.mp4 -i movies/${1}_star_xy.mp4 -i movies/${1}_gas_xy.mp4 -filter_complex "[0:v][1:v]hstack[top];[2:v][3:v]hstack[bottom];[top][bottom]vstack[out]" -map "[out]" movies_stitch/${1}.mp4

ffmpeg -i movies/GSE2iso_fg0.7-lvl4-rng_-15_15_-15_15_BoxCenter_nres1024_gas_xz.mp4 -i movies/GSE4iso_fg0.7-lvl4-rng_-15_15_-15_15_BoxCenter_nres1024_gas_xz.mp4 -i movies/GSE5iso_fg0.7-lvl4-rng_-15_15_-15_15_BoxCenter_nres1024_gas_xz.mp4  -filter_complex hstack=inputs=3 movies_stitch/GSE2_GSE4_GSE5.mp4

ffmpeg -i movies/GSE2iso_fg0.7-lvl4-rng_-15_15_-15_15_BoxCenter_nres1024_gas_xz.mp4 -i movies/GSE2iso_fg0.7-lvl4-fSN-rng_-15_15_-15_15_BoxCenter_nres1024_gas_xz.mp4 -i movies/GSE4iso_fg0.7-lvl4-rng_-15_15_-15_15_BoxCenter_nres1024_gas_xz.mp4 -i movies/GSE5iso_fg0.7-lvl4-rng_-15_15_-15_15_BoxCenter_nres1024_gas_xz.mp4  -filter_complex hstack=inputs=4 movies_stitch/GSE2_fSN_GSE4_GSE5.mp4

ffmpeg -i movies/MW5iso-lvl5-beta05-rng_-80_80_-80_80_Tot_COM_nres512_lowdens_gas_xy.mp4 \
       -i movies/MW5iso-lvl5-beta06-rng_-80_80_-80_80_Tot_COM_nres512_lowdens_gas_xy.mp4 \
       -i movies/MW5iso-lvl5-beta067-rng_-80_80_-80_80_Tot_COM_nres512_lowdens_gas_xy.mp4 \
       -i movies/MW5iso-lvl5-beta07-rng_-80_80_-80_80_Tot_COM_nres512_lowdens_gas_xy.mp4 \
       -i movies/MW5iso-lvl5-beta08-rng_-80_80_-80_80_Tot_COM_nres512_lowdens_gas_xy.mp4 \
       -filter_complex hstack=inputs=5 movies_stitch/MW5iso-lvl5_allbeta_rng_-80_80_-80_80_Tot_COM_nres512_lowdens_gas_xy.mp4

ffmpeg -i movies/MW6iso-lvl5-beta05-dens-rng_-80_80_-80_80_BoxCenter_nres512_lowdens2_gas_xy.mp4 \
       -i movies/MW6iso-lvl5-beta06-dens-rng_-80_80_-80_80_BoxCenter_nres512_lowdens2_gas_xy.mp4 \
       -i movies/MW6iso-lvl5-beta067-dens-rng_-80_80_-80_80_BoxCenter_nres512_lowdens2_gas_xy.mp4 \
       -i movies/MW6iso-lvl5-beta07-dens-rng_-80_80_-80_80_BoxCenter_nres512_lowdens2_gas_xy.mp4 \
       -i movies/MW6iso-lvl5-beta08-dens-rng_-80_80_-80_80_BoxCenter_nres512_lowdens2_gas_xy.mp4 \
       -filter_complex hstack=inputs=5 movies_stitch/MW6iso-lvl5_allbeta-dens-rng_-80_80_-80_80_BoxCenter_nres512_lowdens2_gas_xy.mp4

ffmpeg -i movies/MW6iso-lvl5-beta08-fbar02-dens-rng_-80_80_-80_80_BoxCenter_nres512_lowdens_gas_xy.mp4 \
       -i movies/MW6iso-lvl5-beta08-fbar02-vphi03-dens-rng_-80_80_-80_80_BoxCenter_nres512_lowdens_gas_xy.mp4 \
       -i movies/MW6iso-lvl5-beta08-fbar02-vphi04-dens-rng_-80_80_-80_80_BoxCenter_nres512_lowdens_gas_xy.mp4 \
       -filter_complex hstack=inputs=3 movies_stitch/MW6iso-lvl5-beta08-fbar02_allvphi-dens-rng_-80_80_-80_80_BoxCenter_nres512_lowdens_gas_xy.mp4

ffmpeg -i movies/MW7_GSE4-lvl5-rng_-140_140_-140_140_BoxCenter_nres1024_lowdens_gas_xy.mp4 \
       -i movies/MW7_GSE4-lvl4-rng_-140_140_-140_140_BoxCenter_nres1024_lowdens_gas_xy.mp4 \
       -filter_complex hstack=inputs=2 movies_stitch/MW7_GSE4-lvl5_lvl4-rng_-140_140_-140_140_BoxCenter_nres1024_lowdens_gas_xy.mp4

ffmpeg -i movies/MW7_GSE4-lvl5-rng_-140_140_-140_140_BoxCenter_nres1024_lowdens_gas_xy.mp4 \
       -i movies/MW7_GSE4-lvl5-denscut-rng_-140_140_-140_140_BoxCenter_nres512_lowdens_gas_xy.mp4 \
       -filter_complex hstack=inputs=2 movies_stitch/MW7_GSE4-lvl5-allcuts-rng_-140_140_-140_140_BoxCenter_nres1024_lowdens_gas_xy.mp4

ffmpeg -i movies/MW7_GSE4-lvl5-denscut-Ngb64-rng_-80_80_-80_80_Subhalo_nres512_lowdens_gas_xy.mp4 \
       -i movies/MW7_GSE4-lvl5-denscut-Ngb64-rng_-40_40_-40_40_Subhalo_nres512_lowdens_gas_xy.mp4 \
       -i movies/MW7_GSE4-lvl5-denscut-Ngb64-rng_-30_30_-30_30_Subhalo_nres512_lowdens_gas_xy.mp4 \
       -i movies/MW7_GSE4-lvl5-denscut-Ngb64-steep1-rng_-80_80_-80_80_Subhalo_nres512_lowdens_gas_xy.mp4 \
       -i movies/MW7_GSE4-lvl5-denscut-Ngb64-steep1-rng_-40_40_-40_40_Subhalo_nres512_lowdens_gas_xy.mp4 \
       -i movies/MW7_GSE4-lvl5-denscut-Ngb64-steep1-rng_-30_30_-30_30_Subhalo_nres512_lowdens_gas_xy.mp4 \
       -filter_complex "\
            [0:v][1:v][2:v]hstack=inputs=3[top_row]; \
            [3:v][4:v][5:v]hstack=inputs=3[bottom_row]; \
            [top_row][bottom_row]vstack=inputs=2[v]" \
       -map "[v]" movies_stitch/MW7_GSE4-lvl5-denscut-Ngb64-wsteep-allrng_Subhalo_nres512_lowdens_gas_xy.mp4

