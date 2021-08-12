#!/bin/bash

#
#$ -wd /scratch2/$USER/STAR/BES/build
#$ -N Flow
# -q all.q
#$ -l h=(ncx112|ncx115|ncx117|ncx121|ncx124|ncx12[6-7]|ncx130|ncx132|ncx134|ncx136|ncx138|ncx141|ncx$
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

./build/PIDcomb -i $INPUT_FILE -o $OUTPUT_FILE -m $WORK_MODE -g $ENERGY &>> $LOG_FILE
