#!/usr/bin/perl -w

###############################################################################
#
# Takes the hit / no hit file and checks to see if the 
# no hit entry is actually in hit file.
# If no hit entry is not in hit file, adds it to a new nohits file.
#
# Usage:
# perl hitNohitCheck.pl file.hits file.nohits > file.new_nohits
#
###############################################################################

my $hitFile=$ARGV[0];
my $noHitFile=$ARGV[1];

%hitHash=();
%noHitHash=();
#%intHash=();

open(IN,$hitFile);
while ($lineHit = <IN>) {
	chomp($lineHit);
	my($readTagHit,$strH,$chrH,$locH,$seqH,$confH,$scoreH)=split(/\t/,$lineHit);
	$hitHash{$readTagHit} = $readTagHit;
}
close(IN);

# this file is a fastq, and thus is handled this way
open(IN,$noHitFile);
while($L1 = <IN>) {
	chomp($L1);
        $L2=<IN>;
	$L3=<IN>;
	$L4=<IN>;
	chomp($L2);
	chomp($L3);
	chomp($L4);
	$readTagNoHitProc = substr($L1,1);
	$noHitHash{$readTagNoHitProc}=$readTagNoHitProc;
	if (!exists($hitHash{$readTagNoHitProc})) {
		#print STDERR "$L1\n";
		#print STDERR "$L2\n";
		#print STDERR "$L3\n";
		#print STDERR "$L4\n";
		print "$L1\n";
                print "$L2\n";
                print "$L3\n";
                print "$L4\n";

	}
}
close(IN);

# This must be for debugging, to count a hit was also in the no hit
$count=0;

foreach $key (keys(%hitHash)) {
	if (exists($noHitHash{$key})) {
		$count++;
		#print $key, "\n";
	}
}
#print $count, "\n";

