#!/usr/bin/perl -w

use strict;


my $mut = 1e-8;
my $length = 1e6;
my $N = 1e4;
my $r = 1e-9;
my $theta = 4*$N*$length*$mut;
my $rho = 4*$N*$length*$r;
my $samplesize = 20;
my $pop1 = $samplesize/2;
my $pop2 = $samplesize - $pop1;
my $replications = 1000;
my $mig1prob = 1e-3;
# $mig1prob = 0.7e-3;
my $factor = 0.75;
my $mig1 = $mig1prob * 4 * $N;
my $mig2 = 1.5 * $mig1prob * 4 * $N *$factor;


my $theta1 = $theta;
my $theta2 = $theta * $factor;
my $rho1 = $rho;
my $rho2 = $rho*$factor;

my $cmd1 = "ms $samplesize $replications -t $theta1 -r $rho1 100 -I 2 $pop1 $pop2 0 -em 0 2 1 $mig1 -eM 0.01 0 -ej 1 2 1 > ms1.out";
my $cmd2 = "ms $samplesize $replications -t $theta2 -r $rho2 100 -I 2 $pop1 $pop2 0 -em 0 2 1 $mig2 -eM 0.01 0 -ej 1 2 1 > ms2.out";

system($cmd1);
system($cmd2);


#open(IN2, "ms2.out") or die "Couldn't open ms2.out for input";

sub getavdif{

    my $infile = $_[0];
    
    my @freqs1 = ();
    my @freqs2 = ();
    my $index = 0;
    my $segsites = 0;
    my $repl = 0;
    my $diff1 = 0;
    my $diff2 = 0;
    
    open(IN1, "$infile") or die "Coudln't open ms1.out for input";
    my @res = ();
    $res[$replications - 1] = 0;
    while( defined(my $ln = <IN1>)){
	chomp($ln);
	if($ln =~ /^\s*$/){ next; }
	if($repl == 0 && $ln !~ /\/\//){ next; }
	
	if($ln =~ /\/\//){
	    $repl++;
	    $segsites = 0;
	    @freqs1 = ();
	    @freqs2 = ();
	    $index = 0;
	    next;
	}
	if($ln =~ /segsites:\s*(\d+)/){
	    $segsites = $1;
	    $freqs1[$segsites - 1] = 0;
	    $freqs2[$segsites - 1] = 0;
	    next;
	}
	if($ln =~ /positions/){ next; }
	
	my @v = split(//, $ln);
	for(my $j = 0; $j < @v; ++$j){
	    if($index < $pop1){
		$freqs1[$j] += $v[$j];
	    }else{
		$freqs2[$j] += $v[$j];
	    }
	    
	}
	
	if($index == $samplesize - 1){
	    for(my $j = 0; $j < @freqs1; ++$j){
		$diff1 += abs($freqs1[$j]/$pop1 - $freqs2[$j]/$pop2);	    
		#print $repl, "\t", $freqs1[$j]/$pop1, "\t", $freqs2[$j]/$pop2, "\n";
	    }
	    $res[$repl-1] = $diff1/@freqs1;
	    $diff1 = 0;
	}
	$index++;
    }
    close IN1;
    return(@res);
}

my @dif1 = getavdif("ms1.out");
my @dif2 = getavdif("ms2.out");
my @t = ();
for(my $i =0; $i < @dif1; ++$i){
    print $dif1[$i] - $dif2[$i], "\n";
}
    
    
    
    
    
    
    



    


