#!/usr/bin/perl -w

#	To Run: perl countNTA.pl species files > out
#	SPECIES: hsa / mmu 
# 	J. Baran

$species = shift;
foreach $a (@ARGV) {

   open (FILE,$a) or die;
   $header=<FILE>;
   chomp($header);
   @hparts=split(/\t/,$header);
   $tot=<FILE>;
   $total=0;
   $countout=0;
   $countE=0;
   while (<FILE>) {
      chomp();
      @parts=split(/\t/,$_);
      for ($i=0; $i<scalar(@parts); $i++) {
	 $dat{$hparts[$i]} = $parts[$i];
      }
      $NTA=$dat{Count} -$dat{EM} -$dat{E};
      if ($parts[0] !~ /$species/){
	 $countout+=$NTA;
	 $countE+=$dat{EM};
	 $total+=$dat{Count};
      }
   }
   close(FILE);
   $perc =$countout/$total;
   $pE=$countE/$total;
   print "$a\t$countout\t$pE\t$perc\n";
}
