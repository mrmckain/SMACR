#!/usr/bin/perl -w
use strict;

open my $out, ">", "Read_Mapping_Stats.txt";
print $out "Mismatch\tExperiment\tTotal_Aligned\tTotal_Unaligned\tTotal\n";

my %counts;
my @files = <sRNA_Mapping/*aligned_mismatch*.fsa>;
for my $file (@files){
	$file =~ /(\w+)_(\w+)_mismatch(\d+).fsa/;
	my $exp = $1;
	my $map = $2;
	my $mis = $3;

	open my $tfile, "<", $file;
	while(<$tfile>){
		chomp;
		if(/>/){
			my $count = substr($_,1);
			$counts{$mis}{$exp}{$map}+=$count;
		}
	}
}

for my $mis (sort keys %counts){
	for my $exp (sort {$a cmp $b} keys %{$counts{$mis}}){
		my $total = $counts{$mis}{$exp}{unaligned} + $counts{$mis}{$exp}{aligned};
		print $out "$mis\t$exp\t$counts{$mis}{$exp}{aligned}\t$counts{$mis}{$exp}{unaligned}\t$total\n";
	}
}	
