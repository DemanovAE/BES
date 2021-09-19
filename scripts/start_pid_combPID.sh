#!/bin/bash

INPUT_FILE=$1
WORK_MODE=$2
ENERGY=$3
OUTPUT_DIR=/scratch2/$USER/BES/OUT/${ENERGY}GeV/${WORK_MODE}_${ENERGY}GeV_`date '+%Y%m%d_%H%M%S'`

mkdir -p $OUTPUT_DIR/log
mkdir -p $OUTPUT_DIR/root
mkdir -p $OUTPUT_DIR/sge_error
mkdir -p $OUTPUT_DIR/sge_output

for (( i=0; i < 20; i++ ))
do
  INPUT_FIT=$OUTPUT_DIR/root/$RUNID/First_${ENERGY}GeV_ptBin_$i.root
  OUTPUT_FILE=$OUTPUT_DIR/root/$RUNID/Second_${ENERGY}GeV_ptBin_$i.root
  OUTPUT_LOG=$OUTPUT_DIR/log/$RUNID/ptBin_$i.log
  OUTPUT_O_SGE=$OUTPUT_DIR/sge_output/$RUNID/ptBin_$i.out
  OUTPUT_E_SGE=$OUTPUT_DIR/sge_error/$RUNID/ptBin_$i.err

  qsub -o $OUTPUT_O_SGE -e $OUTPUT_E_SGE run_CombPID_nica.sh $INPUT_FILE $INPUT_FIT $OUTPUT_FILE $ENERGY $WORK_MODE $i 
done