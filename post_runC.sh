#!/bin/bash

#example call: ./post_runC.sh /fullpathToDir/g1Results/

for ARG in "$@"
do
  mkdir $ARG/trash
  DIR=`dirname $ARG`
  BASE=`basename $ARG .results`
  NAME2=CollectLog_$BASE
  RUNDIR=$DIR/$BASE
  bsub -J combine "cat $ARG/3p_summary/* >> $ARG/3p_summary.txt; mv $ARG/3p_summary $ARG/trash"
  bsub -J combine "cat $ARG/5p_summary/* >> $ARG/5p_summary.txt; mv $ARG/5p_summary $ARG/trash"
  bsub -J combine "cat $ARG/shift_summary/* >> $ARG/shift_summary.txt; mv $ARG/shift_summary $ARG/trash"
  bsub -J combine "cat $ARG/lenDist_summary/* >> $ARG/lenDist_summary.txt; mv $ARG/lenDist_summary $ARG/trash"
  bsub -J combine "cat $ARG/ed_summary/* >> $ARG/ed_summary.txt; mv $ARG/ed_summary $ARG/trash"
  bsub -J combine "cat $ARG/*_Shrimp_results.bed >> $ARG/Shrimp_results.bed; mv $ARG/*_Shrimp_results.bed $ARG/trash"
done
