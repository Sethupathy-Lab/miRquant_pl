#!/usr/bin/perl -w
use List::Util qw(min max);
use File::Basename;

my $progdirname = dirname(__FILE__);


$file = $ARGV[0];
if ($ARGV[1] eq 'hsa') {
   $tRNALU=$progdirname . '/resources/hg19_tRNAlu.txt';
   $tRNAbed=$progdirname . '/resources/hg19_tRNA12.bed';
}
else{ 
   $tRNALU=$progdirname . '/resources/mm9_tRNAlu.txt';
   $tRNAbed=$progdirname . '/resources/mm9_tRNA12.bed';
}

open (B12,$tRNAbed) or die;
while(<B12>) 
{
   chomp();
   my($chr,$st,$ed,$nm,$sc,$str,$S,$E,$col,$numB,$bsz,$bst) = split(/\t/,$_);
   $len = $ed -$st;
   if ($numB >1) {

      @bsizes = split(",",$bsz);
      @bstarts = split(",",$bst);
      @{$tRNAinfo{$nm}}[0]=0;
      @{$tRNAinfo{$nm}}[1]=$len;
      @{$tRNAinfo{$nm}}[2]=$bsizes[0];
      @{$tRNAinfo{$nm}}[3]=$bstarts[1];
   }
   else {
      @{$tRNAinfo{$nm}}[0]=0;
      @{$tRNAinfo{$nm}}[1]=$len;
      @{$tRNAinfo{$nm}}[2]=0;
      @{$tRNAinfo{$nm}}[3]=0;

   }
}
close(B12);
open (LU,$tRNALU) or die;
while (<LU>) {
   chomp();
   my($a,$b) = split(/\t/,$_);
   @bnames=split(":",$b);
   foreach $bn(@bnames) {
      $lu{$bn} = $a;
#if (!exists($rlu{$a})) {
# $rlu{$a} = $bn;
#}
#else {
	 push (@{$rlu{$a}},$bn);
#}
   }
}
close(LU);

open(TAB,$file) or die;
$header=<TAB>;
chomp($header);
$total=<TAB>;
chomp($total);
@tparts = split(/\t/,$total);
$mappedCount = $tparts[5];
while(<TAB>) {
  chomp();
  @parts=split(/\t/,$_);
  if ($parts[0] =~ /tRNA/) {
     @name=split("_",$parts[0]);

     $st = $tRNAinfo{$name[0]}[0];
     $ed = $tRNAinfo{$name[0]}[1];
     $ist = $tRNAinfo{$name[0]}[2];
     $ied = $tRNAinfo{$name[0]}[3];

# set up RNA strand ouputs
     $mName="$lu{$name[0]}";
     $info{$mName} = $name[0];
     if (!exists($mx{$mName})) {
	$mx{$mName} =0;
     }
     if (!exists($pu{$mName})) {
	if ($ist ==0) {
	   for ($jj = 0; $jj < $name[3]+3; $jj++){ 
	      $pu{$mName}{$jj} = 0;
	   }
	}
	else {
	   for ($jj=0; $jj< $ist; $jj++) {
	      $pu{$mName}{$jj} =0;
	   }
	   for ($jj=$ist; $jj<= $ied; $jj++) {
	      $pu{$mName}{$jj} =-20;
	   }
	   for ($jj=$ied+1; $jj<$ed+3; $jj++) {
	      $pu{$mName}{$jj} =0;
	   }
	}
     }
     
     if (($name[2]< $ist) || ($name[1] > $ied) || ($ist == 0)) {
	for ($ii = $name[1]; $ii<= $name[2]; $ii++) {
	   if (exists($pu{$mName}{$ii}) && ($pu{$mName}{$ii} != -20)) {
	      $pu{$mName}{$ii} += $parts[5];
	      $mx{$mName} =$pu{$mName}{$ii} if ($pu{$mName}{$ii} > $mx{$mName});
	   }
	}
     }
     $cName="$name[0]";
     $info{$cName} = $name[0];
     if (!exists($mx{$cName})) {
	$mx{$cName} =0;
     }
     if (!exists($sz{$cName})) {
	$sz{$cName} = $name[3];
     }
     if (!exists($pu{$cName})) {
	for ($jj = min($name[1],-40); $jj < 0; $jj++){ 
	   $pu{$cName}{$jj} = -20;
	}
	for ($jj = 0; $jj < $name[3]; $jj++){ 
	   $pu{$cName}{$jj} = 0;
	}
	for ($jj = $sz{$cName}; $jj <= max($name[3]+40,$name[2]); $jj++){ 
	   $pu{$cName}{$jj} = -20;
	}
     }
     for ($ii = $name[1]; $ii <= $name[2]; $ii++) {

	$pu{$cName}{$ii} = 0 if (!exists($pu{$cName}{$ii}) || ($pu{$cName}{$ii}==-20));

	$pu{$cName}{$ii} += $parts[5];
	$mx{$cName} =$pu{$cName}{$ii} if ($pu{$cName}{$ii} > $mx{$cName});
     }
     
  }
}
close(TAB);

open (OUTD, ">jD-PileUp.txt") or die;
open (OUTR, ">jR-PileUp.txt") or die;
foreach $trna (keys(%pu)) {
   @locs = sort {$a <=> $b} keys(%{$pu{$trna}});
   $relCount = 100*($mx{$trna} / $mappedCount);

   $infol = join("\t",@{$tRNAinfo{$info{$trna}}});
   if ($trna =~ /^chr/) {
      print OUTD "$trna\t$relCount\t$infol\t$locs[0]";
   }
   else {
      print OUTR "$trna\t$relCount\t$infol\t$locs[0]";
   }
   foreach $l (@locs) {
      if ($pu{$trna}{$l} != -20){
	 $val =  $pu{$trna}{$l} ;
	 $val = 100* ($pu{$trna}{$l} / $mx{$trna}) if ($mx{$trna}!=0);
	 if ($trna =~ /^chr/) {
	    print OUTD "\t$val";
	 }
	 else {
	    print OUTR "\t$val";
	 }
      }
      else {
	 if ($trna =~ /^chr/) {
	    print OUTD "\t-20";
	 }
	 else{
	    print OUTR "\t-20";
	 }
      }
   }
   if ($trna =~ /^chr/) {
      print OUTD "\n";
   }
   else {
      print OUTR "\n";
   }

}

