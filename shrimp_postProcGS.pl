#!/usr/bin/perl -w

###############################################################################
#
# This script does things.
#
# Usage:
#  perl shrimp_postProcGS.pl GRO-seq_file merge.bed_file SPECIES OUTPUT_FILE? 
use Fcntl qw(:flock SEEK_END SEEK_SET);
use File::Basename;

my $progdirname = dirname(__FILE__);
# read in GS values
$GSLib=shift;                                              # assign GRO-seq varible
$windowFile=shift;                                         # NAME_allGS.bed"
$baseDirName=`dirname $windowFile`;                        # set directory of .bed file
$GENOME=shift;                                             # set species to variable
if ($GENOME eq 'hsa') {                                    # Get necessary table based on genome
   open (MIR, "$progdirname/resources/hsa_table.txt");
}
elsif ($GENOME eq 'mmu') {
   open (MIR, "$progdirname/resources/mmu_table.txt");
}
elsif ($GENOME eq 'cast') {
   open (MIR, "$progdirname/resources/cast_table.txt");
}
else {
   open (MIR, "$progdirname/resources/rno_table.txt");
}
%mirList=();
%mirStrand=();
while($mir=<MIR>) {                                                  # loop through miRs in table
   chomp($mir);                                            
   my ($mName,$mchr,$mSt,$mEd,$mStr,$mSeq,$mHp) = split(/\t/,$mir);  # split line into variables
   $ind = index($mHp, $mSeq);                                        # returns index of mSeq in mHp

   $loc = $mSt + $ind;                                               # set start to actual location of seed
   $loc = $mEd - $ind -length($mSeq) +1 if ($mStr eq "-");           # same as above, but for minus string
   $mchr = uc($mchr);                                                # makes chromosome name uppercase
   $mirList{$mchr}{$loc} = $mName;                                   # puts miR name in chromosome location dict
   $mirStrand{$mchr}{$loc} = $mStr;                                  # puts +/- strand info in chromosome location dict
}

$tagCount=0;
$tagPCount=0;
$zCount=0;
chomp($baseDirName);

foreach $dirL(@ARGV) {                                               # loop through directories in agruments 

   my ($base,$readSize) =split (/_/,$dirL);                          # base = 'readSize', readSize = actual length
   my $dir = $baseDirName . "/" . $dirL;                             # set directory name
   opendir(my $dh, $dir) || die "can't opendir $dir : $!";           # get file list from directory
   @files = grep { /\.out/ && -f "$dir/$_" } readdir($dh);
   closedir $dh;

   %hits=();
   %maps=();
   %tags=();
   $perMatch=(); # Keep track of perfect Matches (to throw away). These are matches to too many loci from bowtie.
   foreach $res (@files) {                                           # for each file in directroy (from above)
      open(RES,"$dir/$res") or die "cant open $res";
      $dummy=<RES>;                                                  # header line
	 while($line = <RES>) {                                      # for line in file
	    chomp($line);   
	    my ($readTag,$window,$strand,$cstart,$cend,$rstart,$rend,$rlen,$score,$estr) = split(/\t/,$line);
	    @misMatches = split( /\d+/,$estr);                       # split estr at digit?
	    if (!exists($maps{$readTag}{$window})) {                 # if window not in dictionary
	       $maps{$readTag}{$window} = $score;                    # set these to values
	       $hits{$window}{$readTag} = $line;
	       $tagCount++;
	    }
	    elsif ($score> $maps{$readTag}{$window}) {               # else if it is in dictionary
	       $maps{$readTag}{$window} = $score;                    # do same, but don't change tagCount
	       $hits{$window}{$readTag} = $line;
	    }
	    unless (scalar(@misMatches) >0) {                        # if editString only #s, no mismatches found
	       $perMatch{$readTag}{$window} = 1;       
	    }
	 }
      close(RES);
      $numTags = scalar(keys(%maps));                                # number of keys in maps dict
      print "\n $res: NumTags 1: $numTags\t";                        # print this info out
   }
   $perMatchCount = 0;
#   foreach $k (keys(%perMatch)) {
#      $perMatchCount ++;
#      foreach $i (keys(%{$maps{$k}})) {
#	 delete ($maps{$k}{$i});
#	 delete ($hits{$i}{$k});
#      }
#      delete ($maps{$k});
#   }
   $numTags2 = scalar(keys(%maps));                                  # These should be the same as above since
   print "NumTags 2: $numTags2 ($perMatchCount)\t";                  # section hashed out
#print "Shrimp Discarded: $perMatchCount, Reads\n"; # Remove perfect Matches
   $tname =`mktemp`;                                       # make temporary variable          
   chomp($tname);
   foreach $k (keys(%maps)) {                              # for all keys in maps
      $maxScore = 0;                                     

      @WindArr = ();
      $maxStr='';
      foreach $i (keys(%{$maps{$k}})) {                    # for all keys in maps[key]
	 my($chr,$loc,$str) = split(/[:\(\)]/,$i);         # get chromosome, location, string info
	 if ($maps{$k}{$i}>$maxScore) {                    # if value is greater than maxscore
	    $maxScore = $maps{$k}{$i};                     # set maxscore to value
	    $maxStr=$str ;                                 # set maxstring to value name
	 }
	 elsif ($maps{$k}{$i} == $maxScore) {              # if not greater than maxscore
	    if ($str =~ /R/) {                             # prefer the RNA strand
	       $maxStr=$str ;                              # set maxstring to value name
	    }
	 }
      }
      #  REMOVE all hits with score less than maxscore. Favor Better matches.
      my @tmpArr1 = keys(%{$maps{$k}});                    # make temporary arrage of maps keys
      foreach $i (keys(%{$maps{$k}})) {                    # loop through maps[key]
	 my($chr,$loc,$str) = split(/[:\(\)]/,$i);         # get chromosome, location, string info
	 if ($maps{$k}{$i} < $maxScore) {                  # if value is less than maxScore
	    delete ($maps{$k}{$i});                        # remove key entry from maps dict
	    delete ($hits{$i}{$k});                        # remove key entry from hits dict
	 }
	 elsif ( ($maps{$k}{$i} == $maxScore) && ($maxStr =~ /R/) && ($str !~ /R/)) {    # if equal to max score and max String not RNA strand and string not RNA strand
	    $logic2=($maxStr=~/R/);
	    $logic3=($str!~/R/);
# delete matches to DNA strand if the max score is on RNA strand
	    delete ($maps{$k}{$i});
	    delete ($hits{$i}{$k});
	 }
	 else {
	    $tmp = $i;
	    $tmp =~ s/\(\+\)/:P/g;
	    $tmp =~ s/\(-\)/:M/g;
	    $tmp =~ s/\(R\+\)/:RP/g;                              # JTB
	    $tmp =~ s/\(R-\)/:RM/g;                               # JTB
	    my($chr,$posA,$posB,$str) = split(/[:-]/,$tmp);       # set chromosome, positions, and orientation
	    my ($rT,$win,$str2,$cst2,$cend2,$rstart2,$rend2,$rlen2,$score2,$estr2) = split(/\t/,$hits{$i}{$k});
	    if (($str eq "P")||($str eq "RP")) {                  # JTB or RP??
	       $str = "+";                                        # plus string
	       $t1 = $posA;                                       # 
	       $posA = $t1 + $cst2;
	       $posB = $t1 + $cend2;
	    }
	    else
	    {
	       $str = "-";
	       $t1 = $posB;
	       $posA = $t1 - $cend2;
	       $posB = $t1 - $cst2;
	    }
	    $chr =~ s/CHR/chr/g;
	    $bedline = join("\t",$chr,$posA,$posB,$i,1,$str);
	    push(@WindArr,$bedline);
	 }
      }
@out = keys(%{$maps{$k}});
# $tname = "shrimpt_" . $readSize;
      $denom = 0;
      if ($GSLib eq "NoGS") {                                        # Ours is currently always NoGS
	 foreach $line (@WindArr) {                                  # for each window
	    my($co,$so,$eo,$wino,$gs,$stro) = split(/\t/,$line);     # GS=1 above
	    $denom += $gs;                                         
	    $tags{$k}{$wino} = $gs;
	 }
      }
      else {                                                         # can skip this
	 $filecont=join("\n",@WindArr);
	 `echo "$filecont">$tname`;
	 @gsRes = `gsBedFilt.pl $GSLib $tname`;
	 foreach $line (@gsRes) {
	    chomp($line);
	    my($co,$so,$eo,$wino,$gs,$stro) = split(/\t/,$line);
	    $denom += $gs;
	    $tags{$k}{$wino} = $gs;
	 }
      }

      foreach $i (keys(%{$maps{$k}})) {                              # for each key in maps[key]
	 if ($denom==0) {                                            # if denominator = 0
	    $denom = 1;                                              # set denom to 1
	    $zCount++;                                               # add 1 to zCount
	 }
	 if (!exists(($tags{$k}{$i})) ) {                            # error check, exits if no tag recorded
	    print "exit: No tag recorded for keys :$k \t$i \t\n";
	    exit;
	 }
	 $tags{$k}{$i} /=$denom;                                     # divided equals
	 $tagPCount+=$tags{$k}{$i};                                  # plus equals
      }
   }                                                                 # END: foreach $k (keys(%maps))
   $numTags = scalar(keys(%maps));                                   # get number of keys in maps dict
   print "NumTags 3: $numTags \n";                                   # prints this info
   `rm $tname`;                                                      # removes tname (hased out earlier)

   $dname = $baseDirName . '/g1Results';                             # maybe makes g1Results
   foreach $k (keys(%hits)) {                                        # for each window in hits dict
      $tmp = $k;                                                     # temp equals window
      $tmp =~ s/\(\+\)/:P/g;                                         # check window orientation
      $tmp =~ s/\(R\+\)/:RP/g;
      $tmp =~ s/\(-\)/:M/g;
      $tmp =~ s/\(R-\)/:RM/g;
      my($chr,$st,$sp,$str) = split(/[:-]/,$tmp);                    # get chromosome, start, stop, string orient
      $dirName = $dname . "/" . $chr;                                # set directory name by chromosome
      `mkdir -p $dirName`;                                           # make chromo folder (& g1Results if not made)
      $fname = $dirName . "/" . $k . ".results";                     # set output file name
      if( length(keys(%{$hits{$k}}))>0){                             # if length of value in hits[key] > 0
	 open (OUT, ">>$fname") or die "cant open $fname";           # open output file name
	 unless (flock(OUT, LOCK_EX | LOCK_NB)) {                    # file lock, do not understand
	    $| = 1;
	    print "Waiting for lock...";
	    flock(OUT, LOCK_EX)  or die "can't lock filename: $!";
	    print "got it.\n"
	 } 
	 seek(OUT, 0, SEEK_END);

	 foreach $h (keys(%{$hits{$k}})) {                           # for each key in hits[key]
	    my ($readTag2,$window2,$strand2,$cstart2,$cend2,$rstart2,$rend2,$rlen2,$score2,$estr2) = split(/\t/,$hits{$k}{$h});
	    $pos = $st + $cstart2;                                   # adjust postion by some metric
	    if ($str =~ /M/) {                                       # if string is minus
	       $pos = $sp - $cend2;                                  # adjust accordingly
	    }
# add a counter here somehow to keep track of the reads that were processed by shrimp2
	    $mir='NA';
	    $winStr = "+";
	    $winStr = "-" if ($str =~ /M/);
	    foreach $p (keys (%{$mirList{$chr}})) {                  # for key in mirList[chromosome]
	       if ($mirStrand{$chr}{$p} eq $winStr) {                # if mirStrand[chromosome position = winstr
		  $d = abs($p - $pos);                               # d = absolute value of p - position
		  if ($d<9){                                         # if d is less than 9
		     $mir = $mirList{$chr}{$p};                      # mir = mir from mirlist
		     last;                                           # equivalent to break in python
		  }
	       }
	    }
	    print OUT "$hits{$k}{$h}$mir\t$tags{$h}{$k}\n";
	 }
	 close(OUT);
      }
   }
}
print "TagCount = $tagCount; zeroCount= $zCount;  pseudoCount=$tagPCount;\n";
# what are these counting???

