#!/usr/bin/perl -w

###############################################################################
#
# Usage: perl normRPMMall.pl PATH_TO_FILE > OUTPUT_NAME
# 
#   PATH_TO_FILE: /small_rna_pipeline/MY_PROJECT/*/TAB_3p_summary.txt
#   OUTPUT_NAME: Any output name
#
###############################################################################

$RUNNTA=0;
$species=shift;
foreach $a (@ARGV) {
   open (FILE,$a) or die;
   $header=<FILE>;
   chomp($header);
   @hparts=split(/\t/,$header);
   $tot=<FILE>;
   $total{$a}=0;
   while (<FILE>) {
      chomp();
      @parts=split(/\t/,$_);
      for ($i=0; $i<scalar(@parts); $i++) {
	 $dat{$hparts[$i]} = $parts[$i];
      }
      $NTA=$dat{Count} -$dat{EM} -$dat{E};
      $total{$a}+=$dat{Count}; # total mapped reads
#	 if ($parts[0] =~ /$species/){ # COMMENT OUT TO DO ALL TYPES OF RES.
	    $datout{$a}{$parts[0]}=$dat{Count};
	    $seedOut{$parts[0]} = $dat{Seed};
	    if ($RUNNTA){
	       $datN{$a}{$parts[0]}=$NTA/$dat{Count};
	       $datT{$a}{$parts[0]}=$dat{T}/$dat{Count};
	       $datA{$a}{$parts[0]}=$dat{A}/$dat{Count};
	    }
	    $list{$parts[0]}=1;
#	 }
   }
   close(FILE);
}
foreach $m (keys(%list)) {
   foreach $ds (keys(%datout)) {
      if (exists($datout{$ds}{$m})) {
	 $re=1000000*($datout{$ds}{$m}/$total{$ds});
	 $datout{$ds}{$m} = $re;
	 if ($re>100) { # CHANGE ME FOR DIFF RPMM THRESH
	    $list2{$m} =1;
	 }
      }
   }
}
print "\t";
foreach $ds (keys(%datout)) {
   print $ds;
   my ($A,$B,$C,$D) = split("/",$ds);
#   print "\t$A";
}
print "\n\t";
foreach $ds (keys(%datout)) {
   my ($A,$B,$C,$D) = split("/",$ds);
#   print "\t$B";
}

if ($RUNNTA) {
   foreach $ds (keys(%datout)) {
      my ($A,$B,$C,$D) = split("/",$ds);
      print "\t$A-NTA";
   }
   foreach $ds (keys(%datout)) {
      my ($A,$B,$C,$D) = split("/",$ds);
      print "\t$A-T";
   }
   foreach $ds (keys(%datout)) {
      my ($A,$B,$C,$D) = split("/",$ds);
      print "\t$A-A";
   }
}
print "\n";
foreach $m (keys(%list2)) {
   print "$m\t$seedOut{$m}";
   foreach $ds (keys(%datout)) {
      $re=0;
      if (exists($datout{$ds}{$m})) {
	 $re=$datout{$ds}{$m};
      }
      print "\t$re";
   }
   if ($RUNNTA) {
      foreach $ds (keys(%datout)) {
	 $NTA=0;
	 if (exists($datout{$ds}{$m})) {
	    $NTA=$datN{$ds}{$m};
	 }
	 print "\t$NTA";
      }
      foreach $ds (keys(%datout)) {
	 $NTAT=0;
	 if (exists($datout{$ds}{$m})) {
	    $NTAT=$datT{$ds}{$m};
	 }
	 print "\t$NTAT";
      }
      foreach $ds (keys(%datout)) {
	 $NTAA=0;
	 if (exists($datout{$ds}{$m})) {
	    $NTAA=$datA{$ds}{$m};
	 }
	 print "\t$NTAA";

      }
   }
   print "\n";
}

