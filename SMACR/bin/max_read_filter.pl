#!/usr/bin/perl -w
use strict;

my %seqs;
$ARGV[0] =~ /(\w+)_trimmed.fq/;
my $exp = $1;
open my $out, ">", "$exp\_reads.cleaned.fq";
open my $file, "<", $ARGV[0];
while(<$file>){
	chomp;
	my $id = $_;
	my $seq = readline($file);
	chomp($seq);
	my $third = readline($file);
	chomp($third);
	my $qual = readline($file);
	chomp($qual);
	if(length($seq) <= 30){
		print $out "$id\n$seq\n$third\n$qual\n";
	}
}
	
