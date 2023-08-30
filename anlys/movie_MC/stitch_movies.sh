
#!/bin/bash

ffmpeg -i movies/${1}_star_xz.mp4 -i movies/${1}_gas_xz.mp4 -i movies/${1}_star_xy.mp4 -i movies/${1}_gas_xy.mp4 -filter_complex "[0:v][1:v]hstack[top];[2:v][3:v]hstack[bottom];[top][bottom]vstack[out]" -map "[out]" movies_stitch/${1}.mp4

