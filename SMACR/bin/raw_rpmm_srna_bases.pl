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
        
        my $prime5 = substr($tarray[4], 0, 1);
	my $lenp = length($tarray[4]);
	my $prime3 = substr($tarray[4], -1);
    	
        $data{$mismatch}{$feature}{prime5}{$prime5}{$tarray[9]}{$tarray[4]}{RPMM}=$tarray[6];
	$data{$mismatch}{$feature}{prime5}{$prime5}{$tarray[9]}{$tarray[4]}{RPM}=$tarray[7];
	$data{$mismatch}{$feature}{prime3}{$prime3}{$tarray[9]}{$tarray[4]}{RPMM}=$tarray[6];
	$data{$mismatch}{$feature}{prime3}{$prime3}{$tarray[9]}{$tarray[4]}{RPM}=$tarray[7];
}
}
my %cleaned;

for my $mis (keys %data){
for my $fea (keys %{$data{$mis}}){
for my $end (keys %{$data{$mis}{$fea}}){
for my $sle (sort keys %{$data{$mis}{$fea}{$end}}){
for my $exp (sort keys %{$data{$mis}{$fea}{$end}{$sle}}){

        #my $count=0;
        
        for my $seq (keys %{$data{$mis}{$fea}{$end}{$sle}{$exp}}){
        #        print "$mis\t$fea\t$end\t$sle\t$exp\t$seq\n";
		$cleaned{$mis}{$fea}{$end}{$sle}{$exp}{RPMM}+=$data{$mis}{$fea}{$end}{$sle}{$exp}{$seq}{RPMM};
		$cleaned{$mis}{$fea}{$end}{$sle}{$exp}{RPM}+=$data{$mis}{$fea}{$end}{$sle}{$exp}{$seq}{RPM};
		$counts{$mis}{$fea}{$end}{$sle}{$exp}++;        
        }
}}}}}
for my $mis (keys %cleaned){
open my $out, ">", "Mismatch$mis\_5prime_Bases_Stats.txt";
open my $out2, ">", "Mismatch$mis\_3prime_Bases_Stats.txt";
print $out "Feature\tExperiment\tBase\tReadsPerMillionMapped\tReadsPerMillions\tCount\n";
print $out2 "Feature\tExperiment\tBase\tReadsPerMillionMapped\tReadsPerMillions\tCount\n";
for my $feature (keys %{$cleaned{$mis}}){
	for my $end (keys %{$cleaned{$mis}{$feature}}){
		for my $sle (sort {$a cmp $b} keys %{$cleaned{$mis}{$feature}{$end}}){
			for my $exp (sort {$a cmp $b} keys %{$cleaned{$mis}{$feature}{$end}{$sle}}){
				if($end eq "prime5"){
					print $out "$feature\t$exp\t$sle\t$cleaned{$mis}{$feature}{$end}{$sle}{$exp}{RPMM}\t$cleaned{$mis}{$feature}{$end}{$sle}{$exp}{RPM}\t$counts{$mis}{$feature}{$end}{$sle}{$exp}\n";
				}
				if($end eq "prime3"){
					print $out2 "$feature\t$exp\t$sle\t$cleaned{$mis}{$feature}{$end}{$sle}{$exp}{RPMM}\t$cleaned{$mis}{$feature}{$end}{$sle}{$exp}{RPM}\t$counts{$mis}{$feature}{$end}{$sle}{$exp}\n"; 
				}
			}
		}
	}
}
}
chdir("../");
