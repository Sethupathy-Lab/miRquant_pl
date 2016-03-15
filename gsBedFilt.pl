#!/usr/bin/perl -w

my $data = shift;
my $keyf = shift;


open(DATA,$data) or die "cant open $data";
my %ref=();
while(my $line = <DATA>) {
   chomp($line);
   my @info = split(/\t/,$line);
#print "$info[0] \t";
   push(@{$ref{$info[0]}},$line);
}
close(DATA);


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
#print "$A - $B : $M\n";
   my $line = $ref{$chr}[$M];
   my @arr = split(/\t/,$line);
   if (($arr[1]<=$loc) && ($arr[2]>=$loc)) {
      if ($arr[3] eq $strand) {
	 $ret = $line;
#	 last;
	 $A=$B+1;
      }
      else {
#print "$arr[3] looking for $strand $M";
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
#print "\n";
if ($ret  eq "NF") {
#print "for:$chr:$loc:$strand Not Found\n";
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
   my $val = $m1 + (($m2-$m1)/$N)*($ref);
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
   if ($val>0){
      print "$chr\t$start\t$Lend\t$en\t$val\t$strand\n" ;
   }
   else {
      $EXmatch ++; 
      print STDERR "$chr\t$start\t$Lend\t$en\t$val\t$strand\n" ;
   }
}

}
