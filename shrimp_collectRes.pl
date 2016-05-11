#!/usr/bin/perl -w
use Fcntl qw(:flock SEEK_END);
use List::Util qw[min max];

$dir = $ARGV[0];                                                           # Directory of results
$bt = $ARGV[0] . ".results";                                               # CHR#.results file
$GENOME = $ARGV[1];                                                        # Species
$outDir=`dirname $bt`;                                                     # Directory for outputs
chomp($outDir);                                                            # Remove new line
$TdiscardCount=0;                                                            
$TCount=0;
$CCount=0;
$libDir= $dir . "/../../";                                                 # Directory containing SAMPLE.stats file
opendir(my $dh, $dir) || die "can't opendir $dir : $!";

# Get list of window result files
@files = grep { /\.results/ && -f "$dir/$_" } readdir($dh);                # Get list of SHRiMP results files
%btWins =();
closedir $dh;
#print "$libDir\n";
opendir(my $dh2, $libDir) || die "can't opendir $libDir : $!";
@filesLib = grep { /LIB.fa/ && -f "$libDir/$_" } readdir($dh2);            # Get windows fasta file
#print "@filesLib\n";
closedir $dh2;
%EM =();
my @bedFile;

# Open bowtie results for this CHR
open (BT,$bt);
while ($line = <BT>) {
   chomp($line);
   my($chrB,$offB,$endB,$windB,$countB,$strB) = split(/\t/,$line);
   $oName = join(":",$offB,$endB);
   $EM{$windB}{$oName} = 0 if (!exists($EM{$windB}{$oName}));
   $EM{$windB}{$oName} +=$countB;

   $winFile = "$windB" . ".results";
   $btWins{$winFile} = 1;
}
foreach $f (@files) {
   $btWins{$f} = 1;
}
close(BT);

# Load mmu_mir annotations for mir locations
if ($GENOME eq 'hsa') {
   open (MIR, "resources/hsa_table.txt");
   $TRNAfile="resources/hg19_tRNA12.bed";
}
elsif ($GENOME eq 'mmu'){
   open (MIR, "resources/mmu_table.txt");
   $TRNAfile="resources/mm9_tRNA12.bed";
}
elsif ($GENOME eq 'cast'){
   open (MIR, "resources/cast_table.txt");
   $TRNAFILE="resources/cast_tRNA12.bed";
}
else {
   open (MIR, "resources/rno_table.txt");
   $TRNAfile="resources/rn4_tRNA12.bed";
}

%mirList=();                                                                             
%mirStrand=();                                                                             
%mirLoc=();
#%mirLen=();
%mirClass=();
while($mir=<MIR>) {
   chomp($mir);
   my ($mName,$mchr,$mSt,$mEd,$mStr,$mSeq,$mHp) = split(/\t/,$mir);
   $ind = index($mHp, $mSeq);
   $loc = $mSt + $ind -1;                                              # MMU Table indexes =+1 of bowtie/shrimp
   $loc = $mEd - $ind -1 if ($mStr eq "-");
   $mchr = uc($mchr);                                                  # Make CHR uppercase
   $mirList{$mchr}{$loc} = $mName;                                     # Store mir by location and chr for later searching
   $mirStrand{$mchr}{$loc} = $mStr;                                    # Store mir strand information
   $mirLoc{$mName} = "$mchr:$loc:$mStr";                               # Contains all location information for miR
   $mirClass{$mName} = 0;                                             
# $mirLen{$mName} = length($mSeq);
}

close(MIR);

#open (TRNA, "$TRNAfile");
#while($trna=<TRNA>) {
#chomp($trna) ;
## Add tRNAs to miR structure for labeling tRNA fragments
#my($tchr,$tst,$ted,$tname,$tscore,$tstr) = split(/\t/,$trna);
#
# $loc = $tst;
# $loc = $ted if ($tstr eq '-');
# $mirList{$tchr}{$loc} = $tname;
# $mirStrand{$tchr}{$loc} = $tstr;
# $mirLoc{$tname} = "$tchr:$loc:$tstr";
# $mirClass{$tname} = 1;
##$mirLen{$tname} = $ted-$tst+1;
#}
#close (TRNA);


# For each window file, process all hits in window
$tfile=`mktemp`;
chomp($tfile);
print "$tfile \n";
foreach $res (keys(%btWins)) {                                         # shrimp results files
   %maps=();                                                           # store edits
   %ExMat =();                                                         # count of exact maps
   %mirs=();                                                           # Name of mir/ potential mir
   %tRNAi = ();                                                        # Name of ovlaping tRNA
   %counts = ();                                                       # counts
   %mirLens = ();                                                      # readlenghts
   %countList = ();
   %eCounts=();                                                        # store edits
   my $tmp = $res;
# Replace +- with PM for splitting
   $tmp =~ s/\(\+\)/:P:/g;                                             # replace - with P
   $tmp =~ s/\(-\)/:M:/g;                                              # replace - with M
   $tmp =~ s/\(R-\)/:RM:/g;                                            # replace - with M for tRNA results
   $tmp =~ s/\(R\+\)/:RP:/g;                                           # replace - with P for tRNA results
# Split parts of name
   my ($CHR, $START, $STOP, $STR,$EXT) = split(/[:-]/,$tmp);           # split edited file name into chromosome, contig start, contig end, string, and extension (.results)
   if ($STR eq "P") {                                                  # next bit is setting sSTR to +, -, R+, & R- if P, M, RP, RM, accordingly
      $sSTR="+";
   }
   elsif($STR eq "RP") {
      $sSTR = "R+";
   }
   elsif($STR eq "M") {
      $sSTR = "-";
   }
   else{
      $sSTR = "R-";
   }

   $filename = "$dir/$res";                                            # sets filename to path_to_file + filename (/proj/seth_lab/smallRNA/PROJECT/SAMPLE/IntermediateFiles/g1Results/CHR/file.results)
#print "$tmp\t$filename\n";
   if (-e $filename) {                                                 # if file exists
      %tRNAlu=();                                                     
      open(RES,"$filename") or die "cant open $filename";              # SHRiMP results file
# Process each read in shrimp result file
      while($line = <RES>) {
	 chomp($line);
	 my ($readTag,$window,$strand,$cstart,$cend,$rstart,$rend,$rlen,$score,$estr,$mir,$pcount) = split(/\t/,$line);  # Split SHRiMP result line
	 my $readkey=join("\t",$window,$strand,$cstart,$cend,$rstart,$rend,$rlen,$score,$estr,$mir);                     # Makes dict key
	 $eCounts{$readkey} = 0 if (!exists($eCounts{$readkey}));      # If not in dictionary, put into dictionary; iterated through below
	 $eCounts{$readkey} +=$pcount;
      }
      close(RES);
      foreach $read (keys(%eCounts)){
	 my ($window,$strand,$cstart,$cend,$rstart,$rend,$rlen,$score,$estr,$mir) = split(/\t/,$read);
	 $pcount=$eCounts{$read};                                      # Set pcount to stored pcount from above (very end of result line)
#print "$read\t$pcount\n";

# Get the genomic coord for out to bed file
	 if (($STR eq "M")||($STR eq "RM")) {                          # Check if minus string
	    $genSTART = $STOP - $cend;                                 # set start coord
	    $genEND = $STOP - $cstart +1;                              # set end coord
	 }
	 else                                                          # If plus string
	 {
	    $genSTART = $START +$cstart-1;                             # set start coord
	    $genEND = $START +$cend ;                                  # set end coord
	 }
	 $chrLC =$CHR;                                                 # Set chromosome
	 $chrLC =~ s/CHR/chr/g;                                        # Make chromosome lowercase
	 $strSYM = "+";                                                # set strSYM to +
	 $strSYM = "-" if (($STR eq "M")||($STR eq "RM"));             # If minus, change strSYM to -
	 $line1 = "$chrLC\t$genSTART\t$genEND\t$estr\t$pcount\t$strSYM\n";     # Set line1 to these values
	 $line = "$chrLC\t$genSTART\t$genEND\t$mir\t$pcount\t$strSYM\n";       # Same as above, but editstring is now mir
	 push (@bedFile,$line1);                                       # Put line1 into bedFile array
	 $lkey="$cstart-$cend-$STR";                                   # Set lkey to chromosomeStart-chromosomeEnd-string
	 if (!exists($tRNAlu{$lkey})){                                 # If lkey is in tRNAlu
	    `echo "$line" > $tfile`;                                   # Add line to output file 'tfile'
	    my @tRNAlist=`windowBed -a $tfile -b $TRNAfile -sm -w 40`; # Get overlap with tRNA, with a 40nt buffer

	    if (scalar(@tRNAlist)>=1){                                 # If there are multiple tRNAs
	       $t = $tRNAlist[0];                                      # Use first tRNA on list
	       my @parts = split(/\t/,$t);                             # 
	       $tlen=$parts[8]-$parts[7];
	       $tst = $parts[1] - $parts[7] ;
	       $tst = $parts[8] - $parts[2] if ($parts[5] eq '-');
	       $ted = $parts[2] - $parts[7]-1;
	       $ted = $parts[8] - $parts[1]-1 if ($parts[5] eq '-');
	       if ($parts[15] > 1) { # number of blocks in bed12
		  @blocksz  = split(",",$parts[16]);
		  @blockst  = split(",",$parts[17]);
		  if (($STR eq "RM")   || ($STR eq "RP")) {
		     $tst += ($blockst[1] - $blocksz[0]) if ($tst > $blocksz[0]);
		     $ted += ($blockst[1] - $blocksz[0]) if ($ted > $blocksz[0]);
		  }
	       }
	       $tRNAlu{$lkey} =  join("_",$parts[9] ,$tst,$ted,$tlen,$sSTR) ; 

#print "$window\t$tRNAlu{$lkey}\t$cstart\t$cend\t$estr\t$tst\t$ted\n";
#print "@parts\n";
	    }
	    else {
	       $tRNAlu{$lkey} = "NA";                                  # No tRNA located
	    }
	 }


	 $loc = $cstart; #1-based
	 $pos = $START + $cstart -1;
	 $pos = $STOP - $cstart +1 if (($STR eq "M")||($STR eq "RM"));  # JTB DIFFERENT FOR minus strand... want 5' locus
	 $len = $cend - $cstart +1;
	 $tName =$tRNAlu{$lkey};
	 if ($tName ne 'NA'){
	    $loc="$cstart:$cend:$sSTR" ;
#print "$res\t$loc\t$tName\n";
	 }
#print "$loc\t$pos\t$estr\n";
	 if ($mir eq "NA") {
	    $mir = join (":",$CHR ,$pos , $STR);
	 }
	 if (!exists($maps{$loc})) {
	    $maps{$loc}= $estr;
	    $mirLens{$loc} = $len;
	    $mirs{$loc} = $mir;
	    $tRNAi{$loc} = $tName;
	    $countList{$loc} =$pcount;
	    $counts{$loc} = $pcount;
# $eCounts{$loc}{$estr} = $pcount;
	 }
	 else {
	    $maps{$loc} = join ("," ,$maps{$loc},$estr);
	    $mirLens{$loc} = join ("," ,$mirLens{$loc},$len);
	    $countList{$loc} =join (",", $countList{$loc},$pcount);
	    $counts{$loc} += $pcount;
#if (!exists($eCounts{$loc}{$estr})){
# $eCounts{$loc}{$estr} = $pcount ;
#}
#else{
# $eCounts{$loc}{$estr} += $pcount ;
#}
	 }
      }
   } # check file exists (Some only have bowtie hits!!!!
   $window ="$CHR:$START-$STOP";
   if ($STR eq "M") {
      $window = $window . "(-)";
   }
   elsif ($STR eq "RM") {
      $window = $window . "(R-)";
   }
   elsif ($STR eq "RP") {
      $window = $window . "(R+)";
   }
   else {
      $window = $window . "(+)";
   }
   #ADD BT hits to lists: 
   foreach $item (keys(%{$EM{$window}})) { # Should be empty for RNA strand
      my ($a,$b) = split(":",$item);

      $a++;  #JTB
      $len = $b - $a +1;  

      
      #print "FIXED BT Index error !!!\t";
# $len = $b - $a if ($STR eq "M");  #JTB CHECK plus strand!!!!

      $pos = $START+$a-1;
      $pos = $STOP - $b if ($STR eq "M");
      $mir = join (":",$CHR ,$pos , $STR);

      $winStr = "+";
      $winStr = "-" if ($STR eq "M");
      $endPos = $pos + $len -1 ;
      $endPos = $pos + $len -1 if($STR eq "M");
      $chrLC =$CHR;
      $chrLC =~ s/CHR/chr/g;
      $line = "$chrLC\t$pos\t$endPos\t$mir\t1\t$winStr\n";
      $line1 = "$chrLC\t$pos\t$endPos\tExact\t$EM{$window}{$item}\t$winStr\n";
      push (@bedFile,$line1);
      `echo "$line" > $tfile`;
      my @tRNAlist=`windowBed -a $tfile -b $TRNAfile -sm -w 40`; #>>> get overlap with tRNA
      $tName ="NA";
      if (scalar(@tRNAlist)>=1){ #use first tRNA on list
	 $t = $tRNAlist[0];
	 my @parts = split(/\t/,$t);
	 $tlen=$parts[8]-$parts[7];
	 $tst = $parts[1] - $parts[7] ;
	 $tst = $parts[8] - $parts[2] if ($parts[5] eq '-');
	 $ted = $parts[2] - $parts[7]-1;
	 $ted = $parts[8] - $parts[1]-1 if ($parts[5] eq '-');
	 $tName = join("_",$parts[9] ,$tst,$ted,$tlen,$winStr) ; # second Name Field
#	       print "E:$window\t$tName\t$a\t$b\t$tst\t$ted\n";
#	       print "E:@parts\n";
#print ">>>>>>>   @parts\n$tName\n";
      }
      foreach $p (keys (%{$mirList{$CHR}})) {    
	 if ($mirStrand{$CHR}{$p} eq $winStr) {
	    $d = abs($p-$pos);
	    if ($d<9) {
	       $mir = $mirList{$CHR}{$p};
	       last;
	    }
	 }
      }
      $a="$a:$b:$winStr" if ($tName ne 'NA');
      if (!exists($maps{$a})) {
	 $maps{$a}="EM";
	 $mirLens{$a} = $len;
	 $mirs{$a} = $mir;
	 $tRNAi{$a} = $tName;
	 $countList{$a} = $EM{$window}{$item};
	 $counts{$a} = $EM{$window}{$item}; # Adding EM from Bowtie to Shrimp COunts
# $eCounts{$a}{EM}=$EM{$window}{$item};
      }
      else {
	 $maps{$a} = join ("," ,$maps{$a},"EM");
	 $mirLens{$a} = join ("," ,$mirLens{$a},$len);
	 $counts{$a} +=$EM{$window}{$item};
	 $countList{$a} =join(",",$countList{$a},$EM{$window}{$item});
      }
   }
#sort hit locations by count values... Process each potential hit from max
#cluster all hits around max and mark as visited.
   @keys = sort {$counts{$b} <=> $counts{$a}} keys( %counts);
   if (scalar(@keys)>0) {
      $myLibFile=$filesLib[0];
      $searchName=$res;
      $searchName =~ s/.results//g;
      @grepResults=`grep -A1 "$searchName" $libDir/$myLibFile`;
      $REFString=$grepResults[1];
      chomp($REFString);
#print ">>>>  $myLibFile $REFString\n";
   }
   %SUMMARY = ();
   %p5 = ();
   %p3 = ();
   %E = ();
   %Shift =();
   %LenDist =();
   %visited = ();
#%countE = ();
   $discardCount=0;

   foreach $kkey (@keys) {

      if ($tRNAi{$kkey} eq 'NA') {
	 $key = $kkey;
      }
      else {
	 ($key,$endKey,$kstr) =split(":",$kkey);
#print "JTB: $tRNAi{$kkey}\t$mirLens{$kkey}\n";
      }

      $TCount += $counts{$kkey};
      if (!exists($visited{$kkey})) {
	 $visited{$kkey} = 1;
	 $offset="NA";
	 $key_loc = $key + $START -1;
	 $key_loc = $STOP - $key if (($STR eq "M")||($STR eq "RM"));
	 $minDistance = "NA";
	 $minDistanceS = "NA";
	 $mNameMir="$CHR:$key_loc:$STR";
	 $mNameMir=$tRNAi{$kkey} unless ($tRNAi{$kkey} eq 'NA');


         foreach $p (keys (%{$mirList{$CHR}})) {    
	    if ($STR eq "P") {
	       $winStr = "+";
	       $dS = $key_loc - $p;
	    }
	    else
	    {
	       $winStr = "-" ;
	       $dS = $p - $key_loc;
	    }
	    if ($mirStrand{$CHR}{$p} eq $winStr) {
	       $d = abs($p-$key_loc);
	       if ($minDistance eq "NA") {
		  $minDistance = $d;
		  $minDistanceS = $dS;
	       }
	       else
	       {
		  if ($d < $minDistance) {
		     $minDistance = $d;
		     $minDistanceS = $dS;
		  }
	       }

	       if ($d<9) {
		  $mNameMir = $mirList{$CHR}{$p};
		  last;
	       }
	    }
	 }


	 $offset = $minDistanceS;
	 if (scalar(keys(%{$mirList{$CHR}}))==0) {
	    $offset = 999; # if there is no miRNAs on the CHR, we don't need to modify the miRNA name to be relative to the offset
	 }
	 if (exists($mirLoc{$mNameMir})){
	    ($mirC,$mirL,$mirS) = split(/:/,$mirLoc{$mNameMir});
	    if ($STR eq "M") {
	       $offset = $mirL - $key_loc
	    }
	    else {
	       $offset = $key_loc - $mirL;
	    }
	 }
	 $mNameMir2=$mNameMir;
	 unless ($offset == 0) { # Getting NA here - Y chrom!!!! JTB!!!
	    $dist = abs($offset);
	    $direction = "+";
	    $direction = "-" if ($offset < 0);
	    if ($dist < 40) {
	       $mNameMir2 = "$mNameMir" . "_$direction" . "_$dist";
	    }
	 }
	 $mName = "$mNameMir2,$tRNAi{$kkey},$offset";
	 @lens = split(/,/,$mirLens{$kkey});
	 @edits = split(/,/,$maps{$kkey});
#	 @edits = keys(%{$eCounts{$kkey}});
#print ">>>! $countList{$kkey}\n";
	 @pCountArray = split(/,/,$countList{$kkey});

	 $numEdits = scalar(@edits);
# $numlens = scalar(@lens);
	 if ($counts{$kkey} > 10) { #min expression for miRNA
	 $CCount+=$counts{$kkey};
	    for ($idx=0; $idx<$numEdits; $idx++) {
	       $e = $edits[$idx];
	       $countE=$pCountArray[$idx];
#print "$key\t$e\t$countE\n";
	       $l = $lens[$idx];
	       $LenDist{$mName}{$l} = 0 if (!exists($LenDist{$mName}{$l}));
	       $LenDist{$mName}{$l}+= $countE;
	       @misMatches = split( /\d+/,$e);
	       $nmisMatches = scalar (@misMatches);
	       if (($e eq "EM")||($nmisMatches==0)) {
		  $ExMat{$mName} = 0 if (!exists($ExMat{$mName}));
		  $ExMat{$mName}+=$countE; # Pcount Val >> pseudocount for Bowtie...done in previous scripts
		  $p5{$mName}{E} = 0 if (!exists($p5{$mName}{E}));
		  $p3{$mName}{E} = 0 if (!exists($p3{$mName}{E}));
		  $E{$mName}{E} =0 if (!exists($E{$mName}{E}));
#print "$e\t$countE\n";
	       }
	       else{
#SPLIT edit string into numbers and letters
		  @mCounts = split(/[NACGT]+/,$e);
		  $nmCounts = scalar (@mCounts);
		  @misMatches = split( /\d+/,$e);
		  $nmisMatches = scalar (@misMatches);
		  $misMatchParse = $nmisMatches;
		  $mmParse=1;

		  if ($nmisMatches >0) { # Exact match case!!!
#print "$e: $nmisMatches\t$misMatches[0]\t$countE\n";
		     if ($misMatches[0] ne "") { # 5p difference
			$p5{$mName}{$misMatches[0]}=0 if (!exists( $p5{$mName}{$misMatches[0]}));
			$p5{$mName}{$misMatches[0]}+=$countE;
			if ($nmisMatches >= $nmCounts) { #3p difference

			   $p3{$mName}{$misMatches[-1]} = 0 if (!exists($p3{$mName}{$misMatches[-1]}));
			   $p3{$mName}{$misMatches[-1]} +=$countE; #Add Pcount val
			      $misMatchParse = $nmisMatches-1;
			}
			else
			{ 
			   $p3{$mName}{E} =0 if (!exists($p3{$mName}{E}));
			   $p3{$mName}{E}+=$countE; #Add Pcount val
			}
			$mmParse=0;
		     }
		     else
		     { 
			$p5{$mName}{E} =0 if (!exists($p5{$mName}{E}));
			$p5{$mName}{E}+=$countE; #Add Pcount val
			if ($nmisMatches > $nmCounts) { #3p difference

			   $p3{$mName}{$misMatches[-1]} = 0 if (!exists($p3{$mName}{$misMatches[-1]}));
			   $p3{$mName}{$misMatches[-1]} +=$countE; #Add Pcount val
			      $misMatchParse = $nmisMatches-1;
			}
			else
			{ 
			   $p3{$mName}{E} =0 if (!exists($p3{$mName}{E}));
			   $p3{$mName}{E}+=$countE; #Add Pcount val
			}
			$mmParse=1;
		     }
		  }
		  else {
		     $p5{$mName}{E} =0 if (!exists($p5{$mName}{E}));
		     $p5{$mName}{E}+=$countE; #Add Pcount val
		  }
#print "$nmisMatches\t$nmCounts\t$misMatchParse\n";
		  for ($i = $mmParse; $i < $misMatchParse; $i++) {
		     if ($mCounts[0] eq "") {
			$p = $i+1;
			for ($j = 1; $j<=$i; $j++) {
			   $p+= $mCounts[$j];
			}
		     }
		     else {
			$p = $i;
			for ($j = 0; $j<$i; $j++) {
			   $p+= $mCounts[$j];
			}
		     }
		     for ($pi=0; $pi< length($misMatches[$i]); $pi++) {
			$refN = substr($REFString,$key + ($p+$pi-1) -1,1);
			$edN=substr($misMatches[$i],$pi,1);
			$eName = join (">",$p+$pi,$refN,$edN);
 print "$e: $kkey: $key,$p,$pi, $eName: $tRNAi{$kkey}\t$REFString\n" if ($tRNAi{$kkey} ne 'NA');
			$E{$mName}{$eName} =0 if (!exists($E{$mName}{$eName}));
			$E{$mName}{$eName}+=$countE; #Add Pcount val
		     }
		  }
		  if ($misMatchParse <= 1) {
		     $E{$mName}{E} =0 if (!exists($E{$mName}{E}));
		     $E{$mName}{E}+=$countE; #Add Pcount val
		  }
	       } #not Exact Match
	       $ExMat{$mName} = 0 if (!exists($ExMat{$mName}));
	       $SUMMARY{$mName} = 0 if (!exists ($SUMMARY{$mName}));
	       $SUMMARY{$mName} +=$countE; #Add Pcount val
		  $l="0";
	       $Shift{$mName}{$l} = 0 if (!exists ($Shift{$mName}{$l}));
	       $Shift{$mName}{$l}+=$countE; #Add Pcount val

	    } #for each edit

	 } #enough hits for primary miRNA
	 else {
#print "Discard $key : $counts{$key}\n";
	    $discardCount += $counts{$kkey} ;
	 }
      } #visited
   } #each key

#print "$res\n";
#PRINT outputs to summary files



###  FROM HERE DOWN IS WRITING OUTPUTS FROM DICTS FOR OUTPUTS
   $TdiscardCount += $discardCount;
   if (! -d $outDir . "/5p_summary") {
       system("mkdir $outDir/5p_summary")
   }
   if (! -d $outDir . "/3p_summary") {
       system("mkdir $outDir/3p_summary")
   }
   if (! -d $outDir . "/ed_summary") {
       system("mkdir $outDir/ed_summary")
   }
   if (! -d $outDir . "/shift_summary") {
       system("mkdir $outDir/shift_summary")
   }
   if (! -d $outDir . "/lenDist_summary") {
       system("mkdir $outDir/lenDist_summary")
   }
   if(1) {
      $fname=$outDir . "/5p_summary/${res}_5p_summary.txt";
      open (P5,">$fname");
      seek (P5, 0, SEEK_END);
      foreach $key (keys (%p5)){
	 print P5 "$key\tCount:$SUMMARY{$key}\tEM:$ExMat{$key}";
	 foreach $edit(keys(%{$p5{$key}})) {
	    $item = join(":",$edit,$p5{$key}{$edit});
	    print P5 "\t$item";
	 }
	 print P5 "\n";
      }
      close (P5);
      $fname=$outDir . "/3p_summary/${res}_3p_summary.txt";
      open (P3,">$fname");
      seek (P3, 0, SEEK_END);
      foreach $key (keys (%p3)){
	 print P3 "$key\tCount:$SUMMARY{$key}\tEM:$ExMat{$key}";
	 foreach $edit(keys(%{$p3{$key}})) {
	    $item = join(":",$edit,$p3{$key}{$edit});
	    print P3 "\t$item";
	 }
	 print P3 "\n";
      }
      close (P3);
      $fname=$outDir . "/ed_summary/${res}_ed_summary.txt";
      open (ED,">$fname");
      seek (ED, 0, SEEK_END);
      foreach $key (keys (%E)){
	 print ED "$key\tCount:$SUMMARY{$key}\tEM:$ExMat{$key}";
	 foreach $edit(keys(%{$E{$key}})) {
	    $item = join(":",$edit,$E{$key}{$edit});
	    print ED "\t$item";
	 }
	 print ED "\n";
      }
      close (ED);
      $fname=$outDir . "/shift_summary/${res}_shift_summary.txt";
      open (SH,">$fname");
      seek (SH, 0, SEEK_END);
      foreach $key (keys (%Shift)){
	 print SH "$key\tCount:$SUMMARY{$key}\tEM:$ExMat{$key}";
	 foreach $edit(keys(%{$Shift{$key}})) {
	    $item = join(":",$edit,$Shift{$key}{$edit});
	    print SH "\t$item";
	 }
	 print SH "\n";
      }
      close (SH);
   $fname=$outDir . "/lenDist_summary/${res}_lenDist_summary.txt";
   open (LD,">$fname");
   seek (LD, 0, SEEK_END);
   foreach $key (keys (%LenDist)){
      print LD "$key\tCount:$SUMMARY{$key}\tEM:$ExMat{$key}";
      foreach $edit(keys(%{$LenDist{$key}})) {
	 $item = join(":",$edit,$LenDist{$key}{$edit});
	 print LD "\t$item";
      }
      print LD "\n";
   }
   close (LD);
   }
}

print "iJTB #<=10_Discards:\t$discardCount\t$TdiscardCount\t$TCount\t$CCount\n";
`rm $tfile`;
$fname="$ARGV[0]_Shrimp_results.bed";
print $fname
open (OUTBED,">>$fname");
seek (OUTBED, 0, SEEK_END);
print OUTBED "@bedFile";
close(OUTBED);
