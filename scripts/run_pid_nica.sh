#!/bin/bash

#
#$ -wd /scratch2/$USER/STAR/BES/build
#$ -N Flow
# -q all.q
#$ -l h=(ncx112|ncx115|ncx117|ncx121|ncx124|ncx12[6-7]|ncx130|ncx132|ncx134|ncx136|ncx138|ncx141|ncx144|ncx150|ncx152|ncx159|ncx16[0-9]|ncx17[0-2]|ncx17[4-6]|ncx18[0-1]|ncx18[4-5]|ncx20[1-3]|ncx20[5-8]|ncx21[1-8]|ncx22[4-8]|ncx23[2-8])
#$ -l h_rt=03:30:00
#$ -l s_rt=03:30:00
#

INPUT_FILE=$1
OUTPUT_FILE=$2
WORK_MODE=$3
ENERGY=$4
LOG_FILE=$5

echo "Arguments:"
echo "Input: $INPUT_FILE"
echo "Output: $OUTPUT_FILE"
echo "Mode: $WORK_MODE"
echo "ENERGY: $ENERGY"
echo "Log: $LOG_FILE"

./FemtoDstAnalyzer_PID -i $INPUT_FILE -o $OUTPUT_FILE -m $WORK_MODE -g $ENERGY &> $LOG_FILE
