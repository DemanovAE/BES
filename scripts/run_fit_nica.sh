#!/bin/bash

#
#$ -wd /scratch2/$USER/STAR/BES/build
#$ -N Fit
# -q all.q
#$ -l h=(ncx112|ncx115|ncx117|ncx121|ncx124|ncx12[6-7]|ncx130|ncx132|ncx134|ncx136|ncx138|ncx141|ncx144|ncx150|ncx152|ncx159|ncx16[0-9]|ncx17[0-2]|ncx17[4-6]|ncx18[0-1]|ncx18[4-5]|ncx20[1-3]|ncx20[5-8]|ncx21[1-8]|ncx22[4-8]|ncx23[2-8])
#$ -l h_rt=03:30:00
#$ -l s_rt=03:30:00
#

INPUT_FILE=$1
INPUT_FIT=$2
OUTPUT_FILE=$3
ENERGY=$4
WORK_MODE=$5
PT_BIN=$6
LOG_FILE=$7

echo "Arguments:"
echo "Input: $INPUT_FILE"
echo "Input histo: $INPUT_FIT"
echo "Output: $OUTPUT_FILE"
echo "Mode: $WORK_MODE"
echo "ENERGY: $ENERGY"
echo "pt bin: $PT_BIN"
echo "Log: $LOG_FILE"

./combPID_Gaus -i $INPUT_FILE -f $INPUT_FIT -o $OUTPUT_FILE -g $ENERGY -m $WORK_MODE -p $PT_BIN &>> $LOG_FILE
