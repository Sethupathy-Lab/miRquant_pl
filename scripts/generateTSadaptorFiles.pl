#!/usr/bin/perl -w
#To use: perl generateTSadaptorFiles.pl Position FULL_PATH_FileList
#       Position: 0-based position in name containing index (File_Cond_AAGGCC_etc.fastq: Position=2)
#	FileList: FULL PATH file list to any unique file in g1Results directory
#	Example: perl generateTSadaptorFiles.pl 2 /proj/seth_lab/smallRNA/PROJ/SAMPLE/*.fastq
# 	J. Baran

$ad_begin='TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC';
$ad_end = 'ATCTCGTATGCCGTCTTCTGCTTG';
$Position=shift;
foreach $a (@ARGV) {
   $dir =`dirname $a`;
   $base = `basename $a`;
   @parts = split(/\./,$base);

   $ext=$parts[$#parts];
   $base = `basename $a $ext`;
   chomp($dir);
   chomp($ext);
   chomp($base);
   chdir $dir or die;
   @parts = split ("_",$base);
   $AD = join ("",$ad_begin, $parts[$Position],$ad_end);

   $file = $dir . "/" . $base . "adaptor";
   print "$file\n";
   open(FILE, ">$file") or die;
   print FILE "$AD\n";
   close(FILE);

}

