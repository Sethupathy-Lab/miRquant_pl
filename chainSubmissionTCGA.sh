#!/bin/bash

#bsub -o logName chainSubmission.sh 24 2 hsa QUAL NoGS /proj/seth_lab/smallRNA/BC_data/Results_1/Jzmicro_1_ATCACG_L008_R1_001.fq

date
# MinRNA (cutof size)
MINrna=14  #change if you want smaller fragments to map to the genome
# Overlap with adapter
OVLAP=$1 #10
NERR=$2  #0
#QThresh=$3
GENOME=$3
QUAL=$4
GSFILE=$5;
ERR_RATE=`echo "$NERR/$OVLAP"|bc -l`
#MAXrna=28
if [[ "$GENOME" == "hsa" ]]
then
   tRNA=resources/hg19_tRNA.bed
   tmRNA=resources/hg19_mature_tRNA_LIB.fa
   BI=$JB_GENOME_DIR/hg19                                                                                                                                   
elif [[ "$GENOME" == "mmu" ]]
then
   tRNA=resources/tRNAs.bed
   tmRNA=resources/mm9_mature_tRNA_LIB.fa
   BI=$JB_GENOME_DIR/mm9                                                                                                                                   
else
   tRNA=resources/rn4_tRNA.bed
   tmRNA=resources/rn4_mature_tRNA_LIB.fa
   BI=$JB_GENOME_DIR/rn4                                                                                                                                   
fi

echo "$GENOME $BI $GSFILE"

IDX=1
shift
shift
shift
shift
shift

for ARG in "$@"
do
  ReadLen=`head -n 2 $ARG | tail -n 1 | wc -c`
#MAXrna=`echo "$ReadLen-1-$OVLAP+$NERR"|bc -l`
  MAXrna=`awk '{if (NR%4==2) { if(length($0)>s){s=length($0)}}} END {print s}' $ARG`
#if [ "$MAXrna" -gt "28" ] 
#  then 
#     echo "Max trimmed size = $MAXrna!  Check your overlap and error settings."
#     exit
#  fi
#MAXrna=50
  EXT=`echo "$ARG" |awk -F . '{print $NF}'`
  DIR=`dirname $ARG`
  BASE=`basename $ARG $EXT`
  DIR_I=$DIR/$BASE/IntermediateFiles
  mkdir -p $DIR_I
  OUT=$DIR_I/${BASE}_O${OVLAP}_E$NERR.fq
  UT=$DIR_I/${BASE}_O${OVLAP}_E${NERR}_UT.fq
  ST=$DIR_I/${BASE}_O${OVLAP}_E${NERR}_ST.fq
  LOGN=Log_${BASE}_O${OVLAP}_E$NERR.txt
  ADFile=$DIR/${BASE}adaptor
#  ADAPTER=`cat $ADFile`
#echo "cutadapt -a $ADAPTER -e $ERR_RATE -O $OVLAP -m $MINrna --untrimmed-output=$UT --too-short-output=$ST -o $OUT $ARG"
#cutadapt -a $ADAPTER -e $ERR_RATE -O $OVLAP -m $MINrna --untrimmed-output=$UT --too-short-output=$ST -o $OUT $ARG

  ln -s $ARG $OUT # replace 

  LIB=$DIR_I/${BASE}_O${OVLAP}_E$NERR
  perl -w preSeperateLib.pl $MINrna $MAXrna $LIB 

  for ((len=$MINrna;len<=$MAXrna;len++))
  do
     A1=${LIB}_$len.fq
     E1=`echo "$A1" |awk -F . '{print $NF}'`
     D1=`dirname $A1`
     jname=b_$len
     B1=`basename $A1 $E1`
     O1=$D1/${B1}hits
     U1=$D1/${B1}noHit
     L1=Log_${B1}_bt.txt

# bowtie for solexa1.3 quality scores (MIN6 resting)
     if [[ $QUAL == 33 ]]
     then
	echo "bowtie -q -nomaqround -a -m 20 --phred33-quals  -n 0 -e 70 -l $len -p 8 --seed=197 --un $U1 $BI $A1 $O1"
	bowtie -q -nomaqround -a -m 20 --phred33-quals    -n 0 -e 70 -l $len -p 1 --seed=197 --un $U1 $BI $A1 $O1
     else
	echo "bowtie -q -nomaqround -a -m 20 --solexa1.3-quals  -n 0 -e 70 -l $len -p 8 --seed=197 --un $U1 $BI $A1 $O1"
	bowtie -q -nomaqround -a -m 20 --solexa1.3-quals    -n 0 -e 70 -l $len -p 1 --seed=197 --un $U1 $BI $A1 $O1
     fi
  done     

  OUTgs=${LIB}_allGS.bed
  rm $OUTgs

  for ((len=$MINrna;len<=$MAXrna;len++))
  do
     A1=${LIB}_$len.hits
     B1=${LIB}_$len.noHit
     B2=${LIB}_new_$len.noHits
     T1=${LIB}_temp_$len.bed

     perl hitNohitCheck.pl $A1 $B1 > $B2
     rm $B1

     awk -v N=$len-1 -F"\t" '{print $3"\t"$4"\t"$4+N"\t"$1" "$5" "$6"\t1\t"$2}' $A1 >>$OUTgs
#awk -v N=$len-1 -F"\t" '{print $3"\t"$4"\t"$4+N"\t"$1" "$5" "$6"\t0\t"$2}' $A1 >$T1
#perl gsBedFiltFQGS.pl ~/mylargefs/PS/min6.groseq.genome.coverage.sort.bed $T1 >> $OUTgs 2>>$B2
 done


   OUTm=${LIB}_merge.bed
   OUTt=${LIB}_tmpJTB.bed
   OUTgse=${LIB}_allGSE.bed
   OUTgseu=${LIB}_allGSEu.bed
   OUTl=${LIB}_LIB.fa

   slopBed -b 0  -i $OUTgs -g $BI.chromSizes >$OUTt
   sort -k1,1 -k2,2n $OUTt  > $OUTgse
   rm $OUTt

   awk -F"\t" '{print $1"\t"$2"\t"$3"\tNAME\t1\t"$6}' $OUTgse | uniq > $OUTgseu
#Add tRNA transcripts to windows
   slopBed -b 40 -i $tRNA -g $BI.chromSizes >> $OUTgseu
   mergeBed -d 65 -s -n -i $OUTgseu > $OUTt
   awk '{print $1"\t"$2"\t"$3"\tNAME\t"$4"\t"$5}' $OUTt >$OUTm                                                                                                                      
   rm $OUTt
   slopBed -b 5  -i $OUTm -g $BI.chromSizes >$OUTt
   mv $OUTt $OUTm

   fastaFromBed -s -fi $BI.fa -bed $OUTm -fo $OUTt 
   tr '[:lower:]' '[:upper:]' < $OUTt > $OUTl
   cat $tmRNA >> $OUTl
   rm $OUTt

   lines=`awk '{x++}END{ print x}' $ARG`
   READS=`echo "$lines/4"|bc -l`
   lines=`awk '{x++}END{ print x}' $OUT`
   CREADS=`echo "$lines/4"|bc -l`
   lines=`awk '{x++}END{ print x}' $UT`
   UREADS=`echo "$lines/4"|bc -l`

#BTH=`awk '{print $4}' ${LIB}_temp_*.bed |sort |uniq |wc -l`
   BTM=`awk '{if (NR%4==1) print $0}' ${LIB}*.noHits |sort |uniq |wc -l`
   GSMatches=`wc -l $OUTgs`
   GSUnique=`awk '{print $4}' $OUTgs |sort|uniq|wc -l`
   WIN=`wc -l $OUTm`
   CHARS=`awk '{ if (NR%2==0) print $0}' $OUTl | wc -m`
#UniqEM=`awk '{ print $4}' ${LIB}_temp_*.bed | sort | uniq | wc -l`
   UniqSentSH=`awk '{if (NR%4==1) print $0}' ${LIB}_new_*.noHits | sort | uniq | wc -l`
   echo "RUN Stats: $ARG"
   echo "         Reads: $READS $ARG"
   echo "     Cut Reads: $CREADS $OUT"
   echo "     UnCut Reads: $UREADS $UT"
   TotCA=`echo "$UREADS+$CREADS" |bc`
   ShortCA=`echo "$READS - $TotCA" |bc`
   echo "________________________________________ Total in Files: $TotCA	Short: $ShortCA"
   echo "    BT Hits: $GSUnique"
   echo "    BT Misses: $BTM"
   TotBt=`echo "$GSUnique+$BTM" |bc`
   echo "________________________________________ Total in Files: $TotBt"
   echo " GS>0  Matches: $GSMatches	unique: $GSUnique"
   echo "________________________________________ "
   echo "       Windows: $WIN"
   echo "         CHARS: $CHARS"

   echo "JBSTATS:\t$ARG\t$CREADS\t$ShortCA\t$GSUnique\t$BTM" 
  for ((len=$MINrna;len<=$MAXrna;len++))
  do
     A1=${LIB}_$len.fq
     rm $A1
     NAME=LOG_$len
     bsub -o $NAME -q idle -n 8 -R "span[hosts=1]" perl -w shrimp_proc.pl $QUAL $len $OUTl ${LIB}_ $GENOME $GSFILE
  done
  bsub -o BTpp.log perl -w bt_postProcSortEMParallel.pl $OUTm $OUTgs $GENOME
   
done     
