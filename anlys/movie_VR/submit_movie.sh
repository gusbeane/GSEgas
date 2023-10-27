#!/bin/bash

export LOGS_DIR=logs-$1-$2-$3
mkdir ${LOGS_DIR}

#SBATCH -o logs/OUTPUT_frames.%j.out
#SBATCH -e logs/ERROR_frames.%j.err

export JOBID=`sbatch -o $LOGS_DIR/OUTPUT.%j.out -e $LOGS_DIR/ERROR.%j.err --array=0-$1 --parsable job_movie-create.sh $2 $3`
sbatch --dependency=afterok:$JOBID job_movie-concat.sh $2 $3

