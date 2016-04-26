#!/usr/bin/perl -w

###############################################################################
#
# Called by chainSubmission.sh wrapper, this script submits the unaligned reads
# from the bowtie alignment (allowing no errors) to SHRIMP for alignment with 
# errors.
#
# Usage:
# perl -w shrimp_proc.pl QUAL SIZE LIB_FILE BASE_NAME GENOME GSFILE
#
###############################################################################

use List::Util qw[min max];
use POSIX;
$Qual=shift;      # Qual offset 33/ 64
$size= shift;     # size of samples
$Lib = shift;     # lib file
$base = shift;    # base name of reads
$GENOME = shift;  # genome used
$GSFILE=shift;    # GRO-seq file (NoGS if not avail)

$reads = $base . "new_$size.noHits";         # output file name

$dirName=`dirname $reads`;                   # directory for output
chomp($dirName);                             
$dname = $dirName . '/readSize_' . $size;    # output directory name
`mkdir $dname`;                              # make output directory
$workingDir=`pwd`;                           # set working directory to present working directory
chomp($workingDir);
chdir $dname or die;                         # change directory to output directory

@baseSeed = (1) x $size;
#print "@baseSeed\n";

# defining the number of mismatches allowed at the 3' end of the seed depending on size of the read
$nzeros =0;
$nzeros =1 if ($size>15);                    # 1 mismatch allowed
$nzeros =2 if ($size>19);                    # 2 mismatch allowed
$nzeros =3 if ($size>23);                    # 3 mismatch allowed
for ($i =0; $i < $nzeros; $i++) {            # Set how many mismatches are allowed based on nzeros
   $r = $i+1;
   $baseSeed[-$r]=0;
}
# What is this doing??
@seedList = ();
for ($i=0; $i< $size-$nzeros; $i++) {
#for ($j=$i+1; $j<$size-$nzeros; $j++){ #7Mar2012: only 1 internal 0
      @seed = @baseSeed;
      $seed[$i] = 0;
# $seed[$j]=0;
      $seedT = join("",@seed);
      push @seedList, $seedT;
#      print "$seedT\n";
#}
}
print "Num Seeds: $#seedList\n";


# Submits jobs to SHRIMP for alignment with mismatches
$shrimpFolder= $ENV{SHRIMP_FOLDER};                     # Sets location of SHRIMP

$nSeeds = scalar(@seedList);                            # Get number of seeds           
$nGroups = ceil($nSeeds/16);                            # Get number of groups (rounded up)
for ($i = 0; $i< $nGroups; $i++) {                      # Loop over number of groups
      $prefix = "group$i";                              # set prefix to group#
      $rseeds = scalar(@seedList);                      # get number of seeds
      $n = min (16, $rseeds);                           # use whatever is smaller, 16 or # of seeds
      $a = 16+$i*16;                                    # 16 + 16 * #
      $b = 1+$i*16;                                     # 1 + 16 * #
      $a = $#seedList+1 if ($a > $#seedList);           # sets a to length of seedlist + 1 if a is > than length of sedlist
      @group =  @seedList[-$a..-$b];
      $seedGroup = join(",", @group);                   # Creates seed group (a seed determines where mismatches can occur when aligning)

      ### Example from SHRiMP readme
      #    $ $SHRIMP_FOLDER/utils/project-db.py --seed 00111111001111111100,\
      #    00111111110011111100,00111111111100111100,00111111111111001100,\
      #    00111111111111110000 --h-flag --shrimp-mode cs hsa.miRNA.cDNA.fa
      #
      #    $ ls hsa.miRNA.cDNA-ls.*    # output from above
      #    hsa.miRNA.cDNA-ls.genome  hsa.miRNA.cDNA-ls.seed.1  hsa.miRNA.cDNA-ls.seed.3
      #    hsa.miRNA.cDNA-ls.seed.0  hsa.miRNA.cDNA-ls.seed.2  hsa.miRNA.cDNA-ls.seed.4

      my $python = $ENV{MYPYTHON} || "";                # set location of python

      ### Split genome to fit ram size ###
      `$python $shrimpFolder/utils/split-db.py --ram-size 46 --prefix $prefix --h-flag --seed $seedGroup $Lib &>$prefix.log`;    # submit shrimp job split-db.py; splits the genome to fit ram size
      $nSplits = `ls -l $prefix*.fa |wc -l`;            # gets the files that were output from above script (the split genome chunks)
      chomp ($nSplits);                                 # remove new line character
      @temp = (12) x scalar(@group);                    # temp = 12 * group list len  (eg 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)((I think))
      $tmpName = join("_",@temp);                       # tempName = temp number list joined by _
      $groupName ="";                                   # set groupName to nothing
      for ($j =1; $j <=$nSplits; $j++) {                # this is logging info 
      $groupName = "$groupName $prefix" . "-46gb-$tmpName" . "seeds-$j" . "of" . $nSplits .  ".fa";
      }
      print "$groupName\n";                             # print groupName

      ### Creates projections of the individual genome pieces, so we avoid having to re-project them every time we start a mapping job (Shrimp README) ###
      `$python $shrimpFolder/utils/project-db.py --shrimp-mode ls --h-flag --seed $seedGroup $groupName &>proj_$prefix.log`;     # submit shrimp job project-db.py

      for ($j =1; $j <=$nSplits; $j++) {                                                           # for each genome chunk
	 $seedName = $prefix . "-46gb-" . $tmpName . "seeds-" .$j . "of" . $nSplits . "-ls";       # set seed name to this long thing
#`$shrimpFolder/bin/gmapper-ls -L $seedName $reads -Q -N 16 -F -q -100 -g -100 -e -10 -f -10 -n 1 --min-avg-qv 25 --qv-offset $Qual --shrimp-format > $prefix.out 2>$prefix.err`;
	 `$shrimpFolder/bin/gmapper-ls -L $seedName $reads -Q -N 16 -F -q -100 -g -100 -e -10 -f -10 -n 1 --qv-offset $Qual --shrimp-format > $prefix.out 2>$prefix.err`;       # submit shrimp job gmapper-ls
                # -L        Load the projection of the i-th chunk of the genome.
                # -Q        The input is fastq format
                # -N 16     Use 16 threads
                # -F        Only process the forward strand of the genome
                # -q -100   The score to open a gap along the read sequence.  Should be negative
                # -g -100   The score to open a gap along the genome sequence. Sould be negative
                # -e -10    The score to extend a gap along the genome sequence.  Should be negative
                # -f -10    The score to extend a gap along the read sequence.  Should be negative
                # -n 1      Create a candidate mapping window for a read based on as little as 1 spaced seed match, aka mapping regions are selected by a single spaced seed match
                # --qv-offset       Interpret qvs in fastq input as PHRED+.  The default is 33 for colour space and 64 for letter space.
                # --shrimp-format   Select old SHRiMP output format.
                #
                # SHRiMP format example
                #     >947_1567_1384_F3       reftig_991      +       22901   22923   3       \
                #         25      25      2020    18x2x3
                #
                #         Additionally, the beginning of  each output file begins  with a specification of
                #         the tab-delimited fields. For example,  the following applies  to the above read
                #         hit:
                #
                #         #FORMAT: readname contigname strand contigstart contigend readstart readend readlength score editstring
                #
                #         The #FORMAT: line allows new fields to be unambiguously added, or others removed or reordered without requiring alteration to the parsing routines.
                #
                #         Descriptions of the columns are as follows:
                #         'readname'  Read tag name
                #         'contigname'    Genome (Contig/Chromosome) name
                #         'strand'    Genome strand ('+' or '-')
                #         'contigstart'   Start of alignment in genome (beginning with 1, not 0).
                #         'contigend' End of alignment in genome (inclusive).
                #         'readstart' Start of alignment in read (beginning with 1, not 0).
                #         'readend'   End of alignment in read (inclusive).
                #         'readlength'    Length of the read in bases/colours.
                #         'score'     Alignment score
                #         'editstring'    Edit string: read to reference summary; see below.
                #         
                #         The edit  string consists  of numbers, characters   and the following additional
                #         symbols: '-', '(' and ')'. It is constructed as follows:
                #         <number> = size of a matching substring
                #         <letter> = mismatch, value is the tag letter
                #         (<letters>) = gap in the reference, value shows the letters in the tag
                #         - = one-base gap in the tag (i.e. insertion in the reference)
                #         x = crossover (inserted between the appropriate two bases)
                #
                #         For example:
                #         A perfect match for 25-bp tags is: "25"
                #         A SNP at the 16th base of the tag is: "15A9"
                #         A four-base insertion in the reference: "3(TGCT)20"
                #         A four-base deletion in the reference: "5----20"
                #         Two sequencing errors: "4x15x6"   (i.e. 25 matches with 2 crossovers)
# Clean up
	 `rm $seedName*seed*`;                          # remove files with this name
      }
      `rm $groupName *.genome`;                         # remove files with this name
}



$merge = $base . 'merge.bed';                           # create name for following submission
`bsub -o PP.$size $workingDir/shrimp_postProcGS.pl $GSFILE $merge $GENOME readSize_$size`
