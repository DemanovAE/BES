#!/bin/bash

#-Set path to the femtoDst files
INPUT_DIR=$1
#-Set the maximum number of files in 1 list
NUM_DIVIDES=$2
#-Set ootput dir

#-Example usage
# . GenerateLists.sh /scratch2/parfenov/StData/27gev/run1/ 100

CURRENT_DIR=${PWD}
mkdir -p ../lists
mkdir -p ../lists/lists27GeV

OUTPUT_DIR="../lists/lists27GeV"

TOTAL_NUM_FILES=`ls $INPUT_DIR/*.femtoDst.root | wc -l`
echo "Total number of DST files: ${TOTAL_NUM_FILES}"

MULT=$((TOTAL_NUM_FILES / NUM_DIVIDES))
RESIDUE=$((TOTAL_NUM_FILES % NUM_DIVIDES))

echo "Mult: $MULT, Residue: $RESIDUE"

for i in `seq 1 $MULT`
do
  OUTPUT_FILE=$OUTPUT_DIR/StRuns${i}.list
  echo $OUTPUT_FILE
  ls $INPUT_DIR/*.femtoDst.root | head -n$((i*NUM_DIVIDES)) | tail -n$NUM_DIVIDES &> $OUTPUT_FILE
done

#-Generate last filelist for residual files
OUTPUT_FILE=$OUTPUT_DIR/StRuns$((i+1)).list
echo $OUTPUT_FILE
ls $INPUT_DIR/*.femtoDst.root |  tail -n$RESIDUE &> $OUTPUT_FILE
