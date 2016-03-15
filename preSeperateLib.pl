#!/usr/bin/perl

###############################################################################
#
# This script is called by the chainSubmission wrapper.  The purpose of this 
# script is to parse the trimmed reads into different size fastqs.
#
# Usage: perl preSeperateLib.pl minRNAlen maxRNAlen cutadapt_outfile_basename
#
###############################################################################

my $short = shift;
my $long = shift;
my $lib_file = shift;

# $long = 28 if ($long > 28);

print "preSeperateLib $short\t$long\t$lib_file\n";

open(LIB,"$lib_file.fq") or die "cant open $lib_file";
my @files;
print "before for\n";
for (my $a = $short; $a<=$long; $a++) {                    # for #'s between minRNAlen and maxRNAlen
   local *FILE;
   open(FILE, ">${lib_file}_$a.fq") or die "foobar ";      # opens a file for reads of some length
   push(@files,*FILE);                                     # adds length file to file list
}
   open(NP, ">${lib_file}_notProc.fq") or die;             # opens a not processed file

print "@files\n";                                          # prints list of files
print "$NP\n";                                             # print not processed file
my $headline = <LIB>;

while($line = <LIB>) {
   chomp($headline);
   chomp($line);
   my $len = length($line);
   my $idx = $len-$short;                                  # gets index of file to append to
   my $file=$files[$idx];                                  # set file to file of that length
#print "$line: $len: $idx\n";
   if (($idx <= $#files)&&($idx>=0)) {                     # make sure index is within read length range

#      print "$headline\n$line\n";
      print $file "$headline\n$line\n";                    # print headline and line to file
      for (my $b=0; $b<2; $b++) {
	 $line2=<LIB>;
	 chomp($line2);
	 print $file "$line2\n";
      }
      $headline=<LIB>;
   }
   else{                                                   # if index not in read length range, add to not proc
      print NP "$headline\n$line\n";
      for (my $b=0; $b<2; $b++) {
	 $line2=<LIB>;
	 chomp($line2);
	 print NP "$line2\n";
      }
      $headline=<LIB>;

   }


}
close (LIB);

foreach my $file (@files)
{
   close $file;
}

