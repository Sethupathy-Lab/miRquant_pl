#!/usr/bin/perl -w

###############################################################################
#
# Usage: perl hist.pl PATH_TO_FQ > OUTPUT_NAME
#
#   PATH_TO_FQ: /sm_rna_pipeline/MY_PROJECT/*/IntermediateFiles/*O10_E1.fq
#   OUTPUT_NAME: Any output name for length distribution
#
###############################################################################

foreach $a (@ARGV) {

open (FILE,$a) or die;
$tot{$a}=0;
while(<FILE>) {
   $seq=<FILE>;
   $sep=<FILE>;
   $qual=<FILE>;

   chomp($seq);
   $len = length($seq);
   $list{$len}=1;
   if (!exists($countI{$a}{$len})) {
      $countI{$a}{$len}=1;
   }
   else {
      $countI{$a}{$len}+=1;
   }
   $tot{$a}++;
}
close(FILE);

}

foreach $a (@ARGV) {
   print "\t$a";
}
print "\n";
foreach $l (keys(%list)) {
   print "$l";
   foreach $a (@ARGV) {
      $re=0;
      if (exists($countI{$a}{$l})) {
	 $re=($countI{$a}{$l}/$tot{$a});
      }
      print "\t$re";
   }
   print "\n";
}
