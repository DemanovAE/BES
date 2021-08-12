#!/bin/bash

INPUT_LIST_DIR=$1
WORK_MODE=$2
ENERGY=$3
OUTPUT_DIR=/scratch2/$USER/BES/OUT/${ENERGY}GeV/CombPID_${WORK_MODE}_${ENERGY}GeV_`date '+%Y%m%d_%H%M%S'`

mkdir -p $OUTPUT_DIR/log
mkdir -p $OUTPUT_DIR/root
mkdir -p $OUTPUT_DIR/sge_error
mkdir -p $OUTPUT_DIR/sge_output

ls -1 $INPUT_LIST_DIR/* | while read line
do
  OUTPUT_ROOT=$OUTPUT_DIR/root/$RUNID/`basename ${line%.*t}`.root
  OUTPUT_LOG=$OUTPUT_DIR/log/$RUNID/`basename ${line%.*t}`.log
  OUTPUT_O_SGE=$OUTPUT_DIR/sge_output/$RUNID/`basename ${line%.*t}`.out
  OUTPUT_E_SGE=$OUTPUT_DIR/sge_error/$RUNID/`basename ${line%.*t}`.err

  qsub -o $OUTPUT_O_SGE -e $OUTPUT_E_SGE start_pid_combPID.sh $line $OUTPUT_ROOT $WORK_MODE $ENERGY $OUTPUT_LOG

done
