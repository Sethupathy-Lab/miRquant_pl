#!/usr/bin/perl -w

###############################################################################
#
# Usage: perl generate_mapping_info_table.pl PATH_TO_FILES
#   PATH_TO_FILES: PATH to .stats files in each unique sample directory
#      Example: /small_rna_pipeline/PROJECT_NAME/*/*.stats
#
# To Run: bsub -o LogName.log perl generate_mapping_info_table.pl PATH_TO_FILES
#
# Output saved as MappingInfoTable.txt
#
# J. Baran-Gale
#
###############################################################################

$startDir=`pwd`;
chomp($startDir);

if(1) {
foreach $a (@ARGV) {
   $dir =`dirname $a`;
   chomp($dir);
   $fdir=join("/",$startDir,$dir);
   chdir $fdir or die("can't move to fdir $dir\n");


   `(head -n 1 TAB_3p_summary.txt && grep -P "chr.*tRNA" TAB_3p_summary.txt | sort -k 6,6nr) > TAB_3p_summary_tRNA.txt`;
   $tRNAcount=0;
   open (FILE,'TAB_3p_summary_tRNA.txt') or die("Can't open file TAB_3p_summary_tRNA.txt\n");
   <FILE>;
   while (<FILE>) {
      chomp();
      @parts = split(/\t/,$_);
      $tRNAcount+=$parts[5];
   }
   close (FILE);

   @dirParts=split(/\//,$dir);
   $Name=pop(@dirParts);
   chop($Name);
  
   `echo "tRNAMapped: $tRNAcount" >> $Name.stats`;
}
}
$outFile=join('/',$startDir,'MappingInfoTable.tsv');
open(OUT,">$outFile");
foreach $a(@ARGV) {
   $dir =`dirname $a`;
   chomp($dir);
   @dirParts=split(/\//,$dir);
   $Name=pop(@dirParts);
   chop($Name);
   $file = join ("/",$startDir,$dir,"$Name.stats");
   open (FILE,$file) or die("Can't open file: $file\n");
   while (<FILE>) {
      chomp();
      my ($A,$B) = split(":",$_);
      $Table{$Name}{$A}=$B;
   }
   close (FILE);
}

print OUT "Sample name";
foreach $N(keys(%Table)){
   print OUT "\t$N";
}
print OUT "\n";

#print Sample Name
$k='file';
print OUT "File name";
foreach $N(keys(%Table)){
    print OUT "\t$Table{$N}{$k}";
} 
print OUT "\n";
#print Total Reads
$k='TotReads';
print OUT "Total reads";
foreach $N(keys(%Table)){
    print OUT "\t$Table{$N}{$k}";
} 
print OUT "\n";
$tr=$k;


#print Trimmed Reads
$k='TrimmReads';
print OUT "Trimmed reads";
foreach $N(keys(%Table)){
    print OUT "\t$Table{$N}{$k}";
} 
print OUT "\n";
#print % trimmed
print OUT "\% Trimmed reads";
foreach $N(keys(%Table)){
    $val=100*$Table{$N}{$k}/$Table{$N}{$tr};
    print OUT "\t$val%";
} 
print OUT "\n";

$tmr=$k;

#print Short Reads
$k='ShortReads';
print OUT "Short reads";
foreach $N(keys(%Table)){
    print OUT "\t$Table{$N}{$k}";
} 
print OUT "\n";
#print % short
print OUT "\% Short";
foreach $N(keys(%Table)){
    $val=100*$Table{$N}{$k}/$Table{$N}{$tr};
    print OUT "\t$val%";
} 
print OUT "\n";


#print Exact Match
$k='EMhits';
print OUT "Exact match to genome";
foreach $N(keys(%Table)){
    print OUT "\t$Table{$N}{$k}";
} 
print OUT "\n";
#print % short
print OUT "\% EM";
foreach $N(keys(%Table)){
    $val=100*$Table{$N}{$k}/$Table{$N}{$tmr};
    print OUT "\t$val%";
} 
print OUT "\n";


#print Non-Exact Match
$k='EMmiss';
print OUT "No exact match to genome";
foreach $N(keys(%Table)){
    print OUT "\t$Table{$N}{$k}";
} 
print OUT "\n";
#print % short
print OUT "\% NEM";
foreach $N(keys(%Table)){
    $val=100*$Table{$N}{$k}/$Table{$N}{$tmr};
    print OUT "\t$val%";
} 
print OUT "\n";


#print Total Mapped
$k='Mapped';
print OUT "Total mapped reads";
foreach $N(keys(%Table)){
    print OUT "\t$Table{$N}{$k}";
} 
print OUT "\n";
#print % 
print OUT "\% Mapped";
foreach $N(keys(%Table)){
    $val=100*$Table{$N}{$k}/$Table{$N}{$tmr};
    print OUT "\t$val%";
} 
print OUT "\n";
$map=$k;


#print Total Mapped to miRs
$k='miRMapped';
print OUT "Total mapped to miRs";
foreach $N(keys(%Table)){
    print OUT "\t$Table{$N}{$k}";
} 
print OUT "\n";
#print % 
print OUT "\% of total mapped to miRs";
foreach $N(keys(%Table)){
    $val=100*$Table{$N}{$k}/$Table{$N}{$map};
    print OUT "\t$val%";
} 
print OUT "\n";


#print Total Mapped to miRs
$k='tRNAMapped';
print OUT "Total mapped to tRNAs";
foreach $N(keys(%Table)){
    print OUT "\t$Table{$N}{$k}";
} 
print OUT "\n";
#print % 
print OUT "\% of total mapped to tRNAs";
foreach $N(keys(%Table)){
    $val=100*$Table{$N}{$k}/$Table{$N}{$map};
    print OUT "\t$val%";
} 
print OUT "\n";


close(OUT);
