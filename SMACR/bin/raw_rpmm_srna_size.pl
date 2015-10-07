#!/usr/bin/env perl -w
use strict;
chdir "Feature_Annotated";
my @files = <*stripped.txt>;

my %counts;
my %mapping;
my %total;


my %data;

for my $tfile (@files){
open my $file, "<", $tfile;
while(<$file>){
        chomp;
        $tfile =~ /FullSet_Mismatch(\d+)_Results_(\w+).txt_stripped.txt/;
        my $mismatch = $1;
        my $feature = $2;

        if(/Chromosome/){
                next;
        }
        
        my @tarray= split/\s+/;
        
        
        $data{$mismatch}{$feature}{$tarray[2]}{$tarray[9]}{$tarray[4]}{"RPMM"}=$tarray[6];
	$data{$mismatch}{$feature}{$tarray[2]}{$tarray[9]}{$tarray[4]}{"RPM"}=$tarray[7];
}
}

my %cleaned;

for my $mismatch (keys %data){
for my $feature (keys %{$data{$mismatch}}){
	for my $size (sort {$a <=> $b} keys %{$data{$mismatch}{$feature}}){
        for my $exp (sort {$a cmp $b} keys %{$data{$mismatch}{$feature}{$size}}){
                for my $seq (keys %{$data{$mismatch}{$feature}{$size}{$exp}}){
                        $cleaned{$mismatch}{$feature}{$size}{$exp}{"RPMM"}+=$data{$mismatch}{$feature}{$size}{$exp}{$seq}{"RPMM"};
                        $cleaned{$mismatch}{$feature}{$size}{$exp}{"RPM"}+=$data{$mismatch}{$feature}{$size}{$exp}{$seq}{"RPM"};
			$counts{$mismatch}{$feature}{$size}{$exp}++;
                }
        }
}
}}

for my $mismatch (keys %cleaned){
open my $out, ">", "Size_ReadCounts_PerExperiment_Mismatch$mismatch.txt";
	print $out "Feature\tExperiment\tSize\tReadsPerMillionMapped\tReadsPerMillion\tTotalCount\n";
for my $feature (keys %{$cleaned{$mismatch}}){
	for my $siz (sort {$a <=> $b} keys %{$cleaned{$mismatch}{$feature}}){
        for my $exp (sort {$a cmp $b} keys %{$cleaned{$mismatch}{$feature}{$siz}}){
                        print $out "$feature\t$exp\t$siz\t$cleaned{$mismatch}{$feature}{$siz}{$exp}{RPMM}\t$cleaned{$mismatch}{$feature}{$siz}{$exp}{RPM}\t$counts{$mismatch}{$feature}{$siz}{$exp}\n";
        }
}
}}
chdir("../");
