#!/usr/bin/perl -w
#To use: perl process_allsummary2tab.pl SCRIPTDIR Species FULL_PATH_FileList
#	SCRIPTDIR: Directory where you perl script is installed (summary2Tab_clust.pl
#	SPECIES: hsa / mmu 
#	FileList: FULL PATH file list to any unique file in g1Results directory
#		Example: /proj/seth_lab/smallRNA/PROJ/SAMPLE/IntermediateFiles/g1Results/shift_summary.txt
# perl process_allsummary2tab.pl /proj/seth_lab/Jeanette/smrnapipeline /proj/seth_lab/smallRNA/PROJ/*/IntermediateFiles/g1Results/shift_summary.txt
# J. Baran
# To Run: bsub -o logfileName.log perl process_allsummary2tab.pl SCRIPTDIR Species FULL_PATH_FileList
$scriptDir=shift;
$species=shift;

print "$scriptDir\n";
$sum2tab= $scriptDir . "/summary2Tab_clust.pl";
foreach $a (@ARGV) {
   $dir =`dirname $a`;
   chomp($dir);
   print "$dir\n";
   chdir $dir or die("can't move to dir $dir\n");

   `perl $sum2tab lenDist_summary.txt 0 $species `;
   `perl $sum2tab 3p_summary.txt 0 $species `;
   `perl $sum2tab ed_summary.txt 0 $species `;

   `(head -n 1 TAB_3p_summary.txt && grep $species TAB_3p_summary.txt | sort -k 6,6nr) > TAB_3p_summary_miR.txt`;
   $CountLine=`head -n 2 TAB_lenDist_summary.txt |tail -n 1`;
   @countParts=split(/\t/,$CountLine);
   $Count = $countParts[5];
   $miRcount=0;
   open (FILE,'TAB_3p_summary_miR.txt') or die("Can't open file TAB_3p_summary_miR.txt\n");
   $head=<FILE>;
   while (<FILE>) {
      chomp();
      @parts = split(/\t/,$_);
      $miRcount+=$parts[5];
   }
   close (FILE);
   `mv TAB_3p_summary_miR.txt TAB_lenDist_summary.txt TAB_3p_summary.txt TAB_ed_summary.txt Shrimp_results.bed ../../`;

   @dirParts=split(/\//,$dir);
   pop(@dirParts);
   pop(@dirParts);
   $topDir=join("/",@dirParts);
   $Name=pop(@dirParts);
   chop($Name);
  
   chdir $topDir or die("can't change to dir $topDir\n");
   `bsub -J tarring tar -zcvf $Name.tgz IntermediateFiles`;
   `echo "Mapped: $Count" >> $Name.stats`;
   `echo "miRMapped: $miRcount" >> $Name.stats`;
}
foreach $a(@ARGV) {
   $dir =`dirname $a`;
   chomp($dir);
   @dirParts=split(/\//,$dir);
   pop(@dirParts);
   pop(@dirParts);
   $topDir=join("/",@dirParts);
   $Name=pop(@dirParts);
   $projDir=join("/",@dirParts);
   chop($Name);
   $file = join ("/",$topDir,"$Name.stats");
   open (FILE,$file) or die("Can't open file: $file\n");
   while (<FILE>) {
      chomp();
      my ($A,$B) = split(":",$_);
      $Table{$Name}{$A}=$B;
      $TKeys{$A}=1;
   }
   close (FILE);
}

foreach $N(keys(%Table)){
   print "\t$N";
}
print "\n";

foreach $k(keys(%TKeys)) {
   print "$k";
   foreach $N(keys(%Table)){
      print "\t$Table{$N}{$k}";
   } 
   print "\n";
}
   

