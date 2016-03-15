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
      $seedGroup = join(",", @group);

      my $python = $ENV{MYPYTHON} || "";                # set location of python
      `$python $shrimpFolder/utils/split-db.py --ram-size 46 --prefix $prefix --h-flag --seed $seedGroup $Lib &>$prefix.log`;    # submit shrimp job split-db.py
      $nSplits = `ls -l $prefix*.fa |wc -l`;            # get number of files ending in .fa (though I don't know why)
      chomp ($nSplits);                                 # remove new line character
      @temp = (12) x scalar(@group);                    # temp = 12 * group list len  (eg 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)((I think))
      $tmpName = join("_",@temp);                       # tempName = temp number list joined by _
      $groupName ="";                                   # set groupName to nothing
      for ($j =1; $j <=$nSplits; $j++) {                # this is logging info 
      $groupName = "$groupName $prefix" . "-46gb-$tmpName" . "seeds-$j" . "of" . $nSplits .  ".fa";
      }
      print "$groupName\n";                             # since groupName is "", not sure why its included
      `$python $shrimpFolder/utils/project-db.py --shrimp-mode ls --h-flag --seed $seedGroup $groupName &>proj_$prefix.log`;     # submit shrimp job project-db.py

      for ($j =1; $j <=$nSplits; $j++) {
	 $seedName = $prefix . "-46gb-" . $tmpName . "seeds-" .$j . "of" . $nSplits . "-ls";       # set seed name to this long thing
#`$shrimpFolder/bin/gmapper-ls -L $seedName $reads -Q -N 16 -F -q -100 -g -100 -e -10 -f -10 -n 1 --min-avg-qv 25 --qv-offset $Qual --shrimp-format > $prefix.out 2>$prefix.err`;
	 `$shrimpFolder/bin/gmapper-ls -L $seedName $reads -Q -N 16 -F -q -100 -g -100 -e -10 -f -10 -n 1 --qv-offset $Qual --shrimp-format > $prefix.out 2>$prefix.err`;       # submit shrimp job gmapper-ls
	 `rm $seedName*seed*`;                          # remove files with this name
      }
#clean up
      `rm $groupName *.genome`;                         # remove files with this name
}



$merge = $base . 'merge.bed';                           # create name for following submission
`bsub -o PP.$size $workingDir/shrimp_postProcGS.pl $GSFILE $merge $GENOME readSize_$size`
