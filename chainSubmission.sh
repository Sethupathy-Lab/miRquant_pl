#!/bin/bash

###############################################################################
#
# Wrapper to run the small RNA seq pipeline.  This is followed by runC.sh
#
# Usage:
#   bsub -o logName chainSubmission.sh 10 1 hsa 33 NoGS /proj/seth_lab/smallRNA/MY_PROJECT/data.fastq
#
###############################################################################

# Import arguments to variables
date
# MinRNA (cutof size)
MINrna=14  #change if you want smaller fragments to map to the genome
# Overlap with adapter
OVLAP=$1 # overlap = 10
NERR=$2  # error = 1
#QThresh=$3
GENOME=$3
QUAL=$4
GSFILE=$5;
ERR_RATE=`echo "$NERR/$OVLAP"|bc -l`
#MAXrna=28

# Check species and retrieve appropriate files
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
elif [[ "$GENOME" == "cast" ]]
then
   tRNA=resources/cast_tRNAs.bed
   tmRNA=resources/cast_mature_tRNA_LIB.fa
   BI=$JB_GENOME_DIR/cast
else
   tRNA=resources/rn4_tRNA.bed
   tmRNA=resources/rn4_mature_tRNA_LIB.fa
   BI=$JB_GENOME_DIR/rn4                                                                                                                                   
fi

echo "Genome = $GENOME Genome path = $BI GS file = $GSFILE"

IDX=1
shift
shift
shift
shift
shift

for ARG in "$@"
do
  ReadLen=`head -n 2 $ARG | tail -n 1 | wc -c`   # gets read length from fastq
  MAXrna=`echo "$ReadLen-1-$OVLAP+$NERR"|bc -l`  # determines the max RNA read length; readLength - 1 - overlap + error
#if [ "$MAXrna" -gt "28" ] 
#  then 
#     echo "Max trimmed size = $MAXrna!  Check your overlap and error settings."
#     exit
#  fi
#MAXrna=50
EXT=`echo "$ARG" |awk -F . '{print $NF}'`        # gets the last field of file name (which is .txt, .csv, ect)
  DIR=`dirname $ARG`                             # directory name of argument
  BASE=`basename $ARG $EXT`                      # combines basename (up to ext) and extension (from above)
  DIR_I=$DIR/$BASE/IntermediateFiles             # points the make directory (next) to correct location with name
  mkdir -p $DIR_I                                # makes intermediatefiles directory
  OUT=$DIR_I/${BASE}_O${OVLAP}_E$NERR.fq         # cutadapt output name
  UT=$DIR_I/${BASE}_O${OVLAP}_E${NERR}_UT.fq     # cutadapt untrimmed output name
  ST=$DIR_I/${BASE}_O${OVLAP}_E${NERR}_ST.fq     # cutadapt too short output name
  LOGN=Log_${BASE}_O${OVLAP}_E$NERR.txt          # Log file name
  ADFile=$DIR/${BASE}adaptor                     # Path to adapter file
  ADAPTER=`cat $ADFile`                          # Adapter sequence
  echo Base name = $BASE
  echo Adaptor file = $ADFile
  echo "cutadapt command = cutadapt -a $ADAPTER -e $ERR_RATE -O $OVLAP -m $MINrna --untrimmed-output=$UT --too-short-output=$ST -o $OUT $ARG"
  cutadapt -a $ADAPTER -e $ERR_RATE -O $OVLAP -m $MINrna --untrimmed-output=$UT --too-short-output=$ST -o $OUT $ARG

#ln -s $ARG $OUT # replace 

  LIB=$DIR_I/${BASE}_O${OVLAP}_E$NERR
  echo LIB variable = $LIB
  echo Running preSeperateLib.pl
  perl -w preSeperateLib.pl $MINrna $MAXrna $LIB 
  echo

#####   BOWTIE ALIGNMENT STEP   #####
# bowtie for solexa1.3 quality scores (MIN6 resting)
#  options: -q = fastq files  -nomaqround = prevents rounding of Phred quality scale  -m 20 = suppress all alignments for read if more than 20 -n 0 = max # of mismatches permitted in seed 
#           -e 70 = max permitted total of quality values at all mismatched read positions -p 8 = permits 8 parallel threads --seed=197 = random number generator --un = unaligned output file

  for ((len=$MINrna;len<=$MAXrna;len++))         # Loop through all possible read lengths
  do
     A1=${LIB}_$len.fq                           # Gets basename from LIB variable and adds read length fastq
     E1=`echo "$A1" |awk -F . '{print $NF}'`     # Gets extension of read length fastq
     D1=`dirname $A1`                            # gets directory name for read length fastq
     jname=b_$len                                # 
     B1=`basename $A1 $E1`                       # gets read length fastq basename (no full path)
     O1=$D1/${B1}hits                            # hits file name
     U1=$D1/${B1}noHit                           # noHit file name
     L1=Log_${B1}_bt.txt                         # Log file which is never used (as far as I can tell)
     echo Read length = $len

     if [[ $QUAL == 33 ]]
     then
	echo "bowtie -q -nomaqround -a -m 20 --phred33-quals  -n 0 -e 70 -l $len -p 8 --seed=197 --un $U1 $BI $A1 $O1"
	bowtie -q -nomaqround -a -m 20 --phred33-quals    -n 0 -e 70 -l $len -p 1 --seed=197 --un $U1 $BI $A1 $O1
     else
	echo "bowtie -q -nomaqround -a -m 20 --solexa1.3-quals  -n 0 -e 70 -l $len -p 8 --seed=197 --un $U1 $BI $A1 $O1"
	bowtie -q -nomaqround -a -m 20 --solexa1.3-quals    -n 0 -e 70 -l $len -p 1 --seed=197 --un $U1 $BI $A1 $O1
     fi
    echo 
  done     

  OUTgs=${LIB}_allGS.bed                         # name for bed file for reads which aligned using bowtie
  rm $OUTgs

  for ((len=$MINrna;len<=$MAXrna;len++))         # Loop through all possible read lengths
  do
     A1=${LIB}_$len.hits                         # Reads which align using bowtie
     B1=${LIB}_$len.noHit                        # Reads which didn't align using bowtie
     B2=${LIB}_new_$len.noHits                   # New no align hits file
     T1=${LIB}_temp_$len.bed                     # Temporary bed file

     perl hitNohitCheck.pl $A1 $B1 > $B2         # double check hit / no hit ??
     rm $B1                                      # Remove bowtie no alignment file

     awk -v N=$len-1 -F"\t" '{print $3"\t"$4"\t"$4+N"\t"$1" "$5" "$6"\t1\t"$2}' $A1 >>$OUTgs      # create bed file from reads which aligned using bowtie
#awk -v N=$len-1 -F"\t" '{print $3"\t"$4"\t"$4+N"\t"$1" "$5" "$6"\t0\t"$2}' $A1 >$T1
#perl gsBedFiltFQGS.pl ~/mylargefs/PS/min6.groseq.genome.coverage.sort.bed $T1 >> $OUTgs 2>>$B2
 done


   OUTm=${LIB}_merge.bed                         # Output name for mergeBed output bed (these are the windows that will be aligned to)
   OUTt=${LIB}_tmpJTB.bed                        # Temp output name (pre-sorted)
   OUTt2=${LIB}_tmpJTB2.bed                       # Temp output name 2
   OUTgse=${LIB}_allGSE.bed                      # Output name for aligned reads bed 
   OUTgseu=${LIB}_allGSEu.bed                    # Output name for aligned reads and tRNAs bed
   OUTl=${LIB}_LIB.fa

# slopBed (from bedtools) increases size of each feature in a bed file be a user-defined number of bases. slopBed restricts resizing to the size of the chromosome (i.e. no start < 0 and no end > chromosome size).
#  options: -b = # of bases to increase by  -i = input file  -g = chomosome sizes

   echo slopBed...
   slopBed -b 0  -i $OUTgs -g $BI.chromSizes >$OUTgse                          # unaltered bed file adjusted for chromosome sizes

   awk -F"\t" '{print $1"\t"$2"\t"$3"\tNAME\t1\t"$6}' $OUTgse | uniq > $OUTt
   sort -k1,1 -k2,2n $OUTt |uniq > $OUTgseu                                    # sorted temp file
   rm $OUTt                                                                    # remove unsorted temp file
#Add tRNA transcripts to windows
   echo Add tRNA transcripts to windows with slopBed...
   slopBed -b 40 -i $tRNA -g $BI.chromSizes >> $OUTgseu                        # add tRNA sequences
   echo Merging files with mergeBed...
   sort -k1,1 -k2,2n $OUTgseu > $OUTt2                                         # sorted temp file
   mergeBed -d 65 -s -c 1 -o count -i $OUTt2 > $OUTt                           # merge bed features that overlap or are withing 65nt of each other
   rm $OUTt2                                                                   # remove temp2 file
   awk '{print $1"\t"$2"\t"$3"\tNAME\t"$5"\t"$4}' $OUTt >$OUTm                 # CHECK THIS IF THINGS GO ASQEW                                                                                                     
   rm $OUTt                                                                    # remove original (non-awked) mergeBed file
   slopBed -b 5  -i $OUTm -g $BI.chromSizes >$OUTt                             # add 5nt to each side of merged bed file (this is the windows to align to)
   mv $OUTt $OUTm                                                              # change temp name to _merge.bed name

   echo fastaFromBed...
   fastaFromBed -s -fi $BI.fa -bed $OUTm -fo $OUTt                             # converts _merge.bed file to fasta to SHRIMP align too (the windows fasta)
   tr '[:lower:]' '[:upper:]' < $OUTt > $OUTl                                  # converts all lowercase to uppercase
   cat $tmRNA >> $OUTl                                                         # append tRNA sequences to windows fasta
   rm $OUTt                                                                    # remove temp

   lines=`awk '{x++}END{ print x}' $ARG`                                       # gets total lines from fastq file
   READS=`echo "$lines/4"|bc -l`                                               # calculates reads from lines (1 read per 4 lines in fastq)
   lines=`awk '{x++}END{ print x}' $OUT`                                       # gets total lines from cutadapt trimmed fastq file
   CREADS=`echo "$lines/4"|bc -l`                                              # calculates trimmed reads from lines
   lines=`awk '{x++}END{ print x}' $UT`                                        # gets total lines from cutadapt untrimmed fastq file
   UREADS=`echo "$lines/4"|bc -l`                                              # calculates untrimmed reads from lines

#BTH=`awk '{print $4}' ${LIB}_temp_*.bed |sort |uniq |wc -l`
   BTM=`awk '{if (NR%4==1) print $0}' ${LIB}*.noHits |sort |uniq |wc -l`       # number of unaligned reads
   GSMatches=`wc -l $OUTgs`                                                    # number of aligned reads
   GSUnique=`awk '{print $4}' $OUTgs |sort|uniq|wc -l`                         # number of unique aligned reads
   WIN=`wc -l $OUTm`                                                           # number of windows
   CHARS=`awk '{ if (NR%2==0) print $0}' $OUTl | wc -m`                        # number of characters in window fasta
#UniqEM=`awk '{ print $4}' ${LIB}_temp_*.bed | sort | uniq | wc -l`
   UniqSentSH=`awk '{if (NR%4==1) print $0}' ${LIB}_new_*.noHits | sort | uniq | wc -l`            # Not sure on this
   echo "RUN Stats: $ARG"
   echo "         Reads: $READS $ARG"
   echo "     Cut Reads: $CREADS $OUT"
   echo "     UnCut Reads: $UREADS $UT"
   TotCA=`echo "$UREADS+$CREADS" |bc`                                          # trimmed reads + untrimmed reads
   ShortCA=`echo "$READS - $TotCA" |bc`                                        # total reads - (trimmed reads + untrimmed reads)
   echo "________________________________________ Total in Files: $TotCA	Short: $ShortCA"
   echo "    BT Hits: $GSUnique"
   echo "    BT Misses: $BTM"
   TotBt=`echo "$GSUnique+$BTM" |bc`
   echo "________________________________________ Total in Files: $TotBt"
   echo " GS>0  Matches: $GSMatches	unique: $GSUnique"
   echo "________________________________________ "
   echo "       Windows: $WIN"
   echo "         CHARS: $CHARS"

   StatsFile=$DIR/$BASE/${BASE}stats             # Create stats file
   echo "file:$ARG">$StatsFile                   # stats file full path
   echo "TotReads:$READS">>$StatsFile            # total reads
   echo "TrimmReads:$CREADS">>$StatsFile         # trimmed reads
   echo "ShortReads:$ShortCA">>$StatsFile        # short reads
   echo "EMhits:$GSUnique" >> $StatsFile         # EMhits
   echo "EMmiss:$BTM">>$StatsFile                # EMmisses
  for ((len=$MINrna;len<=$MAXrna;len++))         # loops over all lengths
  do
     A1=${LIB}_$len.fq
     rm $A1                                      # remove length fastq
     NAME=LOG_$len                               # set NAME to LOG_length file
#bsub -x -o $NAME perl -w shrimp_proc.pl $QUAL $len $OUTl ${LIB}_ $GENOME $GSFILE
    bsub -o $NAME -q idle -n 8 -R "span[hosts=1]" perl -w shrimp_proc.pl $QUAL $len $OUTl ${LIB}_ $GENOME $GSFILE
 #    bsub -x -o $NAME                                perl -w shrimp_proc.pl $QUAL $len $OUTl ${LIB}_ $GENOME $GSFILE
  done
  echo "bsub -o BTpp.log perl -w bt_postProcSortEMParallel.pl $OUTm $OUTgs $GENOME"
  bsub -o BTpp.log perl -w bt_postProcSortEMParallel.pl $OUTm $OUTgs $GENOME
   
done     
