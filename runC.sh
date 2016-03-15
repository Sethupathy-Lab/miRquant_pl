#!/bin/bash

#example call: ./runC.sh hsa /fullpathToDir/g1Results/CHR*.results

GENOME=$1
shift

for ARG in "$@"
do
  DIR=`dirname $ARG`
  BASE=`basename $ARG .results`
  NAME2=CollectLog_$BASE
  RUNDIR=$DIR/$BASE
  bsub -o $NAME2 -q idle perl -w shrimp_collectRes.pl $RUNDIR $GENOME
#  bsub -o $NAME2 perl -w shrimp_collectRes.pl $RUNDIR $GENOME
done
