#!/usr/bin/perl -w

###############################################################################
#
# Usage: perl genNormalRPM.pl SPECIES PATH_TO_FILES > OUTPUT_NAME 
#   SPECIES: hsa / mmu
#   PATH_TO_FILES: /small_rna_pipeline/PROJECT_NAME/*/TAB_3p_summary.txt
#   OUTPUT_NAME: Any output name you choose
#
# J. Baran
###############################################################################

$species = shift;
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
#$NTA=$dat{Count} -$dat{EM} -$dat{E};
      $total{$a}+=$dat{Count}; # total mapped reads
#      if ($parts[0] =~ /$species/){     # this is for getting only miRs, leave hashed for all
	 $datout{$a}{$parts[0]}=$dat{Count};
	 $list{$parts[0]}=1;
#      }
   }
   close(FILE);
}

foreach $ds (keys(%datout)) {
   print "\t$ds";
}
   print "\n";
foreach $m (keys(%list)) {
   print "$m";
   foreach $ds (keys(%datout)) {
      $re=0;
      if (exists($datout{$ds}{$m})) {
	 $re=1000000*($datout{$ds}{$m}/$total{$ds});
      }
      print "\t$re";
   }
   print "\n";
}

