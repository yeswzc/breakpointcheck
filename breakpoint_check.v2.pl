#!/usr/bin/perl -w
#
#Author: Zhichao Wu
#
use strict;
use warnings;

my $vcf = shift or die "Usage: $0 <vcf> <bam>\n";
my $bam = shift or die "Usage: $0 <vcf> <bam>\n";

my $minMQ=40;
#average mapping depth
my ($avgDP,$n)=(0,0);
open DP,"samtools depth -Q $minMQ $bam |" || die $!;
while(<DP>){
	my @e=split;
	$n++;
	$avgDP+=$e[2];
	last if $n>=10000;
}
close DP;
$avgDP=sprintf "%.2f",$avgDP/$n;
print "###AverageHighQualMapDepth:$avgDP\n";



###map insert size and variation==================================================================================
#=================================================================================================================
my %insCK;
my %ins;
$n=0;
open SM,"samtools view $bam -q $minMQ |" || die $!;
while(<SM>){
	next if /^@/;
	$n++;
	my @e=split;
	next unless $e[1] & 66;
	next if $e[5] =~ /S/;
	next if $e[4] < $minMQ;
	next if $e[8] <= 100;
	my $rg="unset";
	for(@e[11..$#e]){
		if(/^RG/){
			$rg=(split /:/,$_)[-1];
			last;
		}
	}
	push @{$insCK{$rg}},$e[8];
	my @tmp;
	next if $n<20000;
	foreach (keys %insCK){
		push @tmp,scalar @{$insCK{$_}};
	}
	my @ttmp=grep $_>10000,@tmp;
	last if scalar @ttmp==scalar @tmp;
}
close SM;

foreach my $rg (keys %insCK){
	my @ins=sort {$a<=>$b} @{$insCK{$rg}};
	@ins=@ins[5000..$#ins-5000];
	my $sum=0;
	$sum+=$_ for @ins;
	my $mean=sprintf "%.2f", $sum/scalar(@ins);
	my $sqsum = 0;
	$sqsum+=( ($_-$mean)**2 ) for @ins;
	my $sd =sprintf "%.2f", sqrt($sqsum/scalar(@ins));
	$ins{$rg}{'mean'}=$mean;
	#$ins{$rg}{'upper'}=sprintf "%.2f",$mean+3*$sd;
	#$ins{$rg}{'lower'}=sprintf "%.2f",$mean-3*$sd;
	$ins{$rg}{'upper'}=sprintf "%.2f",$mean*1.5;
	$ins{$rg}{'lower'}=sprintf "%.2f",$mean/1.5;
	print "###RG:$rg; Mean:$mean; SD:$sd; Upper:$ins{$rg}{'upper'}; Lower:$ins{$rg}{'lower'};\n";
}
%insCK=();

###read VCF and Check=============================================================================================
#=================================================================================================================
if($vcf=~/.gz$/){
	open IN,"gzip -dc $vcf |" or die $!;
}else{
	open IN, $vcf or die $!;
}
while (my $line=<IN>) {
	if ($line=~/^#/){print $line;next;}
	chomp $line;
	my @e = split /\s+/, $line;
	my $fh;
	my $chr1 = $e[0];
	my $end = 0;
	if ($line =~ /END=(\d+);/) {
		$end = $1;
	}
	my $chr2=$chr1;
	if ($line =~ /CHR2=(\S+?);/) {
		$chr2 = $1;
	}
	#bp reads covrage
	my $covbp1=bp_dep($chr1,$e[1]);
#next if ($covbp1>5*$avgDP);
	my $covbp2=bp_dep($chr2,$end);
#next if ($covbp2>5*$avgDP);
	my $cipos=($1+$2)/2 if($line=~/CIPOS=-?(\d+),(\d+)/); #for Delly
	$cipos=1000 if $cipos<1000;
	#check reads around bp
	#1. multiple mate chromosomes?
	#2. many peecentage of reads with low mapping quality?
	#3. High variation of map insert size?
	my $beg = $e[1]-$cipos; my $fin = $e[1]+$cipos;
	$fin=$beg+2*$cipos if $fin<0;
	my $cmd1="samtools view $bam $chr1:$beg-$fin";
	my ($lowQualRate,$upper,$lower,$readN,%mate)=bp_check($cmd1);
	print "$line\tBP1;COV:$covbp1;LowQualRate:$lowQualRate;INSlower:$lower;INSupper:$upper;READN:$readN";
	my @mt=grep $mate{$_}>=2,keys %mate;
	if (@mt){
		print "\tBP1MATE";
		print ";$_:$mate{$_}" for @mt;
	}
	@mt=();
	($lowQualRate,$upper,$lower,$readN,%mate)=();
	my $beg2 = $end-$cipos; my $fin2 = $end+$cipos;
	$fin2=$beg2+2*$cipos if $fin2<0;
	my $cmd2="samtools view  $bam $chr2:$beg2-$fin2";
	($lowQualRate,$upper,$lower,$readN,%mate)=bp_check($cmd2);
	print "\tBP2;COV:$covbp2;LowQualRate:$lowQualRate;INSlower:$lower;INSupper:$upper;READN:$readN";
	@mt=grep $mate{$_}>=3,keys %mate;
	if (@mt){
		print "\tBP2MATE";
		print ";$_:$mate{$_}" for @mt;
	}
	print "\n";
}
close IN;

1;
#==============================================================================================
#==============================================================================================

#reads coverage, no matter maping quality
sub bp_dep{
	my ($chr,$bp)=@_;
	my ($a,$b)=($bp-100,$bp+100);
	$a=0 if $a<0;
	my ($dp,$n)=(0,0);
	my $cmd="samtools depth $bam -r $chr:$a-$b";
	open DP,"$cmd |" || die $!;
	while(<DP>){
		my @e=split;
		$n++;
		$dp+=$e[2];
	}
	close DP;
	$n=1 if $n==0;
	$dp=sprintf "%.2f", $dp/$n;
	return $dp;
}

#===============================================================================================
#===============================================================================================
#check reads around bp
#1. multiple mate chromosomes?
#2. many peecentage of reads with low mapping quality?
#3. High variation of map insert size?
sub bp_check{
	my $cmd=shift;
	my ($readN,$highQualReadN,$lowQual)=(0,0,0);
	my ($upperN,$lowerN)=(0,0);
	my %mate;
	open SM,"$cmd |" || die $!;
	while(<SM>){
		next if /^@/;
		my @e=split;
		$readN++;
		$lowQual++,next if $e[4]<10;
		$highQualReadN++;
	
		my $rg;
		for(@e[11..$#e]){
			if(/^RG/){
				$rg=(split /:/,$_)[-1];
				last;
			}
		}
		
		if( $e[6] ne "="){
			next if $ins{$rg}{'mean'}>=1500; #this script is used for ins=~460bp, in rice 3k RG
			$mate{$e[6]}++;
		}
		next if $e[8]<100; ###
		next unless defined $ins{$rg};
		next if $ins{$rg}{'mean'}>=1500; #this script is used for ins=~460bp, in rice 3k RG
		my $upper=$ins{$rg}{'upper'};
		my $lower=$ins{$rg}{'lower'};
		$upperN++ if($e[8])>=$upper;
		$lowerN++ if($e[8])<=$lower;
	}
	close SM;
	$readN=1 unless $readN;
	my $lowQualRate=sprintf "%.3f",$lowQual/$readN;
	return($lowQualRate,$upperN,$lowerN,$highQualReadN,%mate);
}

sub median {
	my @vals = @_;
	my $len = @vals;
	if ($len % 2) {
		return $vals[int($len/2)];
	} else {
		return ($vals[$len/2] + $vals[$len/2-1])/2;
	}
}

sub sum {
	my @vals = @_;
	my $ret = 0;
	foreach (@vals) {
		$ret += $_;
	}
	return $ret;
}


