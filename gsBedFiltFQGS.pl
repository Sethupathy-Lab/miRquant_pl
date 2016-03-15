#!/usr/bin/perl -w

my $data = shift;
my $keyf = shift;


my %ref=();
open(DATA,$data) or die "cant open $data";
while(my $line = <DATA>) {
   chomp($line);
   my @info = split(/\t/,$line);
#print "$info[0] \t";
   push(@{$ref{$info[0]}},$line);
}
close(DATA);


%Visited =();
%arrMatch =();
%arrNoHit = ();
$EXmatch=0;
open (KEYF,$keyf) or die;
while ($key=<KEYF>) {
   chomp($key);
   my ($chr,$start,$Lend,$en,$score,$strand) = split(/\t/,$key);
   my $loc = $start + int(($Lend-$start)/2);


   my $A = 0;
   my $B = $#{$ref{$chr}};
#print "$A - $B \t";

   my $ret = "NF";
   while ($A <=$B) {
      my $M = int(($A+$B)/2);
      my $line = $ref{$chr}[$M];
      my @arr = split(/\t/,$line);
      if (($arr[1]<=$loc) && ($arr[2]>=$loc)) {
	 if ($arr[3] eq $strand) {
	    $ret = $line;
	    $A=$B+1;
	 }
	 else {
	    $A = ($strand eq "-") ? $A-1: $A+1;
	    $B = ($strand eq "-") ? $B-1: $B+1;
	 }
      }
      elsif ($loc < $arr[1]) {
	 $B = $M-1;
      }
      else {
	 $A = $M+1;
      }
   }
   if ($ret  eq "NF") {
# IF not Found: (usually chrM or chrY) Put back in Hit anyway!
      my ($tag,$seq,$qual) = split(/\s/,$en);
      print "$chr\t$start\t$Lend\t$tag\t1\t$strand\n" ; # hits file in pipeMIN6.sh
	 $EXmatch ++; 
      $arrMatch{$tag} = 1;
#if (!exists($arrNoHit{$tag})) {
# $headLn = "@" . $tag;
#@ASCII = unpack("C*",$qual);
#for ($i = 0; $i<=$#ASCII; $i++) {
# $ASCII[$i] += (64-33);
#}
# $qual = pack("C*",@ASCII);
#
# $Visited{$tag} = join("\n",$headLn,$seq,"+",$qual);;
# $arrNoHit{$tag} = 1;
#}
   }
   else {
      my $N = 2000; # defined in dataset
	 my ($ch,$st,$end,$str,$d1,$d2,$As,$Aa,$Is,$Ia) = split(/\t/,$ret);
#my $nw=$end-$st+1;
      my $offset = $loc - $st;
      my $iSkip = int($offset / $N)+1; # 104 val for 102 points
	 my $a = $st + ($iSkip-1)*$N;

      my @dat = split(/[:-]/,$Is);


      if ($strand eq "+") {
	 if ($loc > $a +($N/2)) {
	    $m1 = $dat[$iSkip];
	    $m2 = $dat[$iSkip+1];
	    $d = -$N/2;
	 }
	 else{
	    $m1 = $dat[$iSkip-1];
	    $m2 = $dat[$iSkip];
	    $d = $N/2;
	 }
	 $ref = $loc - $a +$d;
      }
      else {
# - strand runs from end to st.
	 $offset = $end - $loc;
	 $iSkip = int($offset / $N) +1;  # 104 val for 102 points
	    $b = $end - ($iSkip-1)*$N;
	 if ($loc > $b -($N/2)) {
	    $m1 = $dat[$iSkip-1];
	    $m2 = $dat[$iSkip];
	    $d = $N/2;
	 }
	 else{
	    $m1 = $dat[$iSkip];
	    $m2 = $dat[$iSkip+1];
	    $d = -$N/2;
	 }
	 $ref = $b-$loc + $d;
	 $a=$b;
      }
      my $val = int($m1 + (($m2-$m1)/$N)*($ref));
      my $dsz = @dat;
# $t1 = $m2-$m1;
# $t2 = $t1/$N;
# $t3 = $ref;

#print "KJL: $val =  $m1  -  $m2, st: $a, p:$loc\t diff=$t1, m=$t2 , $t3\n ";

#my $val = $dat[$iSkip];
#print "$#dat \t $iSkip \t $val\n";
#print "@dat\n";
#print "for: $chr : $loc : $strand - $ch:$st - $end $str = $val\n";
# $k = $nw/2000;

#print "$nw, $k - $dsz\n";
#print "for:$chr:$loc:$strand  = $val\n";
      my ($tag,$seq,$qual) = split(/\s/,$en);

# last modified 12 Apr 2012 (KAL)
      if ($val>0){
	 print "$chr\t$start\t$Lend\t$tag\t$val\t$strand\n" ; # hits file in pipeMIN6.sh
	    $EXmatch ++; 
	 $arrMatch{$tag} = 1;
      }
      else { # noHits files in pipeMIN6.sh
#if (!exists($Visited{$tag})) 
	 if (!exists($arrNoHit{$tag})) {
	    $headLn = "@" . $tag;
	    @ASCII = unpack("C*",$qual);
	    for ($i = 0; $i<=$#ASCII; $i++) {
	       $ASCII[$i] += (64-33);
	    }
	    $qual = pack("C*",@ASCII);

	    $Visited{$tag} = join("\n",$headLn,$seq,"+",$qual);;
	    $arrNoHit{$tag} = 1;
	 }

#print STDERR "$chr\t$start\t$Lend\t$en\t$val\t$strand\n" ;
      }

# add in lines that will compare @arrMatch and %arrNoHit and delete those in %arrNoHit and %Visited that are in %arrMatch
   }

}
foreach $tag (keys(%arrNoHit)) {
   if (!exists ($arrMatch{$tag})) {
      print STDERR "$Visited{$tag}\n";
   }
}
# change file that is referenced in pipeMIN6.sh so that this is incorporated
