#!/usr/bin/perl -w
use strict;
`mkdir Experiment_Filtered_Results`;
chdir "Experiment_Filtered_Results";
my @files = <../Reformatted_Initial_Results/*sRNA_Results.txt>;

for my $file (@files){
	open my $tfile, "<", $file;
	my $prefix;
	$file =~ /Results\/(.+)\.txt/;
	$prefix = $1;
	open my $outfile, ">", $prefix . "_ExpFiltered.txt";
	my $header = readline($tfile);
	print $outfile "$header";
	my %results;
	my %sub_results;
	
	while(<$tfile>){
		chomp;
		my @tarray = split /\s+/;
		$results{$tarray[0]}{$tarray[1]}{$tarray[2]}{$tarray[9]}= "$tarray[3] $tarray[4] $tarray[5] $tarray[6] $tarray[7] $tarray[8]";
		
		$tarray[9] =~ /(.+)(\d?)/;
		my $exp_prefix = $1;
		my $rep = $2;
		unless($rep){
			$rep = 0;
		}
		$sub_results{$tarray[0]}{$tarray[1]}{$tarray[2]}{$exp_prefix}{$2}=1;
	}

	for my $chr (keys %results){
		for my $start (sort {$a<=>$b} keys %{$results{$chr}}){
			for my $size (sort {$a<=>$b} keys %{$results{$chr}{$start}}){
				for my $experiment (sort keys %{$results{$chr}{$start}{$size}}){
					$experiment =~ /(.+)\d?/;
					if (scalar keys %{$sub_results{$chr}{$start}{$size}{$1}} > 0){
						my @temp = split(" ", $results{$chr}{$start}{$size}{$experiment});
						print $outfile "$chr\t$start\t$size\t$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$experiment\n";
					}
					else{
						next;
					}
				}
			}
		}
	}
}
		
chdir "../";	
