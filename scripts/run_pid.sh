#!/bin/bash

#$ -N Flow

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

./build/FemtoDstAnalyzer_PID -i $INPUT_FILE -o $OUTPUT_FILE -m $WORK_MODE -g $ENERGY &> $LOG_FILE
