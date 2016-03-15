#!/usr/bin/perl -w
use File::Basename;
my $progdirname = dirname(__FILE__);

$res = $ARGV[0];
$ORDER = $ARGV[1];
$GENOME=$ENV{JB_GENOME_DIR};

if ($ARGV[2] eq 'hsa') {
   $mmuFile=$progdirname . '/resources/hsa_tableL.bed';
   $tRNAFile=$progdirname . '/resources/hg19_tRNA.bed';
   $CODE ='hsa';
   $genome =$GENOME . '/hg19.fa';
   $refAnn=$progdirname . '/resources/hg19_ref.bed';
#	$refInt=$progdirname . '/resources/hg19_introns.bed';
}
elsif ($ARGV[2] eq 'mmu') {
   $mmuFile=$progdirname . '/resources/mmu_tableL.bed';
   $tRNAFile=$progdirname . '/resources/tRNAs.bed';
   $CODE ='mmu';
   $refAnn=$progdirname . '/resources/mm9_ref.bed';
#	$refInt=$progdirname . '/resources/mm9_introns.bed';
   $genome =$GENOME . '/mm9.fa';
}
else {
   $mmuFile=$progdirname . '/resources/rno_tableL.bed';
   $tRNAFile=$progdirname . '/resources/rn4_tRNA.bed';
   $CODE ='rno';
   $refAnn=$progdirname . '/resources/rn4_ref.bed';
#	$refInt=$progdirname . '/resources/rn4_introns.bed';
   $genome =$GENOME . '/rn4.fa';
}
$ALL=1;
$ALL = $ARGV[3] if (scalar(@ARGV)==4);

open(MMU,$mmuFile) or die "cant open $mmuFile\n";;
while (<MMU>) {
   chomp();
   my($chr,$st,$ed,$name,$score,$strand) = split(/\t/,$_);

   my($mir,$seq) = split(/:/,$name);
   $bedInfo{$mir} = join("\t",$chr,$st,$ed,$strand);
}
close(MMU);
open(TRNA,$tRNAFile) or die "cant open $tRNAFile\n";
while (<TRNA>) {
   chomp();
   my($chr,$st,$ed,$name,$score,$strand) = split(/\t/,$_);

   $bedInfo{$name} = join("\t",$chr,$st,$ed,$strand);
}
close(TRNA);
%expression =();
%baseExp =();
%baseKeys =();
%bedLine =();
%counts = ();
%features =();
%Novelarr =();
my $tmp = $res;
$tfile=`mktemp`;
chomp($tfile);
open (TBED,">$tfile");
open(RES,"$res") or die "cant open $res";
while($line = <RES>) {
   chomp($line);
   my @info = split (/\t/,$line);
   my $N = $info[0];
   my ($en,$et,$off) = split(/,/,$N);
   my $Nbase = $en;
   $Nbase =~ s/_[+-]_\d//g;
   my ($C,$nC) = split(/:/,$info[1]);
   $expression{$N} = $nC;
   if ($N=~ m/$CODE/) { #if ($Nbase =~ /hsa/) 
      $baseExp{$Nbase} = 0 if (!exists($baseExp{$Nbase}));
      $baseExp{$Nbase} +=$nC;
      push (@{$baseKeys{$Nbase}},$N);
      my($c,$s,$e,$r) = split(/\t/,$bedInfo{$Nbase});
      if ($r eq "+") {
	 $bedLine{$N} = join("\t",$c,$s+$off,$s+$off+7,$N,1,$r);
      }
      else {
	 $bedLine{$N} = join("\t",$c,$e-$off-8,$e-$off-1,$N,1,$r);
      }
   }
   elsif ($N=~ m /tRNA/){
      @parts = split("_",$N);
      $tNbase = $parts[0];
      $o1 = $parts[1];
#JTB chr16.tRNA24-GlyGCC_40_56_71_-
      $baseExp{$tNbase} = 0 if (!exists($baseExp{$tNbase}));
      $baseExp{$tNbase} +=$nC;
      push (@{$baseKeys{$tNbase}},$N);
      my($c,$s,$e,$r) = split(/\t/,$bedInfo{$tNbase});
      if ($r eq "+") {
	 $bedLine{$N} = join("\t",$c,$s+$o1+1,$s+$o1+8,$N,1,$r);
      }
      else {
	 $bedLine{$N} = join("\t",$c,$e-$o1-8,$e-$o1-1,$N,1,$r);
      }
   }
   else 
   {
      my($c,$s,$r) = split(/:/,$Nbase);
      $c =~ s/CHR/chr/g;
      $c =~ s/UN/Un/g;
      $c =~ s/RANDOM/random/g;
      $c =~ s/.chr/chr/g;
      $R = substr($r,0,1);
      if ($R eq "P") {
	 $bedLine{$N} = join("\t",$c,$s+1,$s+8,$N,1,"+");
      }
      else {
	 $bedLine{$N} = join("\t",$c,$s-8,$s-1,$N,1,"-");
      }
      $key="$c:$R";
#push(@{$Novelarr{$key}},$N);
      $Novelarr{$key}{$N}=1;;
      
# $baseExp{$NbaseNov} +=$nC;
#push (@{$baseKeys{$NbaseNov}},$N);
   }
   print TBED "$bedLine{$N}\n";
   for ($i=1; $i<=$#info; $i++) {
      my ($L,$V) = split(/:/,$info[$i]);
      $features{$N}{$L} = $V;
      $tot{$L} = 0 if (!exists ($tot{$L}));
      $tot{$L} += $features{$N}{$L};
      if (!exists($counts{$L})) {
	 $counts{$L} = $V ;
      }
      else {
	 $counts{$L} += $V ;
      }
   }
}
close(RES);
foreach $k (keys(%Novelarr)) {
   @hits = sort { $expression{$b} <=> $expression{$a} } keys(%{$Novelarr{$k}});
   %NovelMax=();
   foreach $h (@hits) {
      my ($en,$et,$off) = split(/,/,$h);
      my $Nbase = $en;
      $Nbase =~ s/_[+-]_\d//g;
      my($c,$s,$r) = split(/:/,$Nbase);
      $R = substr($r,0,1);
      $NbaseNov="$c:$s:$R";
# $num=scalar(keys(%NovelMax));
#print "$num\t";
#print "$NbaseNov\t $s\t ";
      $nkeys = scalar(keys(%NovelMax));
      if ($nkeys<1) {
	 my @keylist=keys(%NovelMax);
#print "first: @keylist\t";
	 $NovelMax{$s} = 1;
      }
      else {
	 $found=0;
	 foreach $item (keys(%NovelMax)) {
	    $distance = abs($s - $item);
	    if ($distance < 9) {
	       $NbaseNov="$c:$item:$R";
#print "$h: $NbaseNov\n";
#print "$item\t";
	       $found=1;
	       last;
	    }
	 }
	 if ($found==0) {
	    $NovelMax{$s} = 1;
#print "unfound\t";
	 }
      }
      $baseExp{$NbaseNov} = 0 if (!exists($baseExp{$NbaseNov}));
      $baseExp{$NbaseNov} +=$expression{$h};
#print " $baseExp{$NbaseNov}\n";
      push (@{$baseKeys{$NbaseNov}},$h);
   }
}
close (TBED);
$tfile2=`mktemp`;
chomp($tfile2);
print "$tfile\n";
`fastaFromBed -fi $genome -s -name -tab -bed $tfile -fo $tfile2`;
#`rm $tfile`;
open(TBED2,$tfile2) or die;
while (<TBED2>) {
   my($N,$seq) = split(/\t/,$_);
   chomp($seq);
   $seed{$N} = uc($seq);
}
close(TBED2);
#`rm $tfile2`;
@ann=`windowBed -w 0 -sm -a $tfile -b $refAnn`;
print "windowBed -w 0 -sm -a $tfile -b $refAnn\n";
foreach $A (@ann) {
   chomp($A);
   @parts = split(/\t/,$A);
   $gAnn{$parts[3]} = $parts[9];
}
@AS=`windowBed -w 0 -Sm -a $tfile -b $refAnn`;
foreach $S (@AS) {
   chomp($S);
   @parts = split(/\t/,$S);
   $gAS{$parts[3]} = $parts[9];
}

if ($ORDER ==1) { #NUM
   @keys = sort { $a <=> $b } keys(%counts);
}
elsif ($ORDER ==2) { #alpha
   @keys = sort { lc($a) cmp lc($b) } keys(%counts);
}
else { # $ORDER ==0 #byVal
   @keys = sort { $counts{$b} <=> $counts{$a} } keys(%counts);
   print "order = 0\n";
}
my @expKeys = sort { $expression{$b} <=> $expression{$a}} keys (%expression);
my @expBaseKeys = sort { $baseExp{$b} <=> $baseExp{$a}} keys (%baseExp);
$fOUT = "TAB_" .$res;
open (OUT,">$fOUT");
print OUT "Name\ttRNA\tmiRbaseOffset\tSeed\tPercent";
foreach $k (@keys) {
   print OUT "\t$k";
}
print OUT "\n";
print OUT "Total\t\t\t\t";
#foreach $e (@expKeys) {
#foreach $k (@keys) {
# $tot{$k} = 0 if (!exists ($tot{$k}));
# $tot{$k} += $features{$e}{$k};
#}
#}
foreach $k (@keys) {
   print OUT "\t$tot{$k}";
}
print OUT "\n";
foreach $e1 (@expBaseKeys) {
#print "$e1\n";
   my @expKeys1 = sort {$expression{$b} <=> $expression{$a}} @{$baseKeys{$e1}};
   my $sumE=0;
   foreach $e (@expKeys1) {
      $sumE+=$features{$e}{Count};
   }
   foreach $e (@expKeys1) {
      my ($en,$et,$off) = split(/,/,$e);
      my $frac = $features{$e}{Count}/$sumE;
      if ($et eq 'NA') {
	 if (exists($gAnn{$e})){
	    $et = $gAnn{$e} ;
	 }
	 else {
	    $et = "AS:" . $gAS{$e} if (exists($gAS{$e}));
	 }
      }
      if (($ALL == 1) || ($en !~ m/CHR/)) {
	 if (!exists($seed{$e})) {
	    $seed{$e} = '';
	 }
	 print OUT "$en\t$et\t$off\t$seed{$e}\t$frac";

	 foreach $k (@keys) {
	    $features{$e}{$k}=0 if (!exists($features{$e}{$k}));
	    print OUT "\t$features{$e}{$k}";
	 }
	 print OUT "\n";
      }
      else {
#	 print "ALL\t$ALL\t$en\n";
      }
   }
}
print "$keys[0]: $tot{$keys[0]}\t$keys[1]: $tot{$keys[1]}\n";

close(OUT);
