#!/usr/bin/perl -w
use strict;
`mkdir Reformatted_Initial_Results`;
chdir "Reformatted_Initial_Results";
my @files = <../sRNA_Mapping/*Results.txt>;

my %data;
my %read_counts;
my $tm;
#my $tm = localtime();
#print "$tm\n";
for my $tfile (@files){
	open my $file, "<", $tfile;
	readline($file);
	$tm = localtime();
	print "$tm: $tfile\n";
	my ($raw_count, $aligned_reads);
	while(<$file>){
		$tm = localtime();
		#print "$tm\n";
		my @tarray = split /\t/;
		$tfile =~ /\_(Mismatch\d)\_/;
		my $mismatch = $1;
		$mismatch =~ /Mismatch(\d)/;
		my $mismatch_count = $1;
		$tfile =~ /sRNA_Mapping\/(\w+)_sRNA_Dmel6.0.1/;
		my $exp = $1;
		unless (exists $read_counts{$exp}{$mismatch}){
			my ($raw_count, $aligned_reads);
			open my $reduced_reads, "<", "../sRNA_Mapping/$exp\.srnas.fsa";
			while(<$reduced_reads>){
				chomp;
				if(/^>/){
					$raw_count += substr($_, 1);
				}
			}
			open my $hit_reads, "<", "../sRNA_Mapping/$exp\_aligned_mismatch$mismatch_count.fsa";
			while(<$hit_reads>){
				chomp;
				if(/^>/){
					$aligned_reads += substr($_, 1);
				}
			}	 
			$read_counts{$exp}{"raw"}=$raw_count;
			$read_counts{$exp}{$mismatch}=$aligned_reads;
		}
		$data{$mismatch}{$tarray[0]}{$tarray[1]}{$tarray[2]}{$exp}= "$tarray[3] $tarray[4] $tarray[5] $tarray[6] $tarray[7]";
	}
}
$tm = localtime();
print "$tm\n";
for my $mismatch (sort keys %data){
	open my $out, ">", "FullSet_" . $mismatch . "_sRNA_Results.txt";
	for my $region (sort keys %{$data{$mismatch}}){
		print $out "Chromosome\tStart\tSize\tDirection\tSeq\tHits\tReadsPerMillionMapped\tReadsPerMillion\tRPMPerHits\tExperiment\n";
		
			for my $start (sort {$a<=>$b} keys %{$data{$mismatch}{$region}}){
				for my $size (sort {$a<=>$b} keys %{$data{$mismatch}{$region}{$start}}){
					for my $exp (sort keys %{$data{$mismatch}{$region}{$start}{$size}}){
					my @narray = split(" ", $data{$mismatch}{$region}{$start}{$size}{$exp});
					my $rpm_mapped = $narray[3]*$read_counts{$exp}{"raw"}/$read_counts{$exp}{$mismatch};
					print $out "$region\t$start\t$size\t$narray[0]\t$narray[1]\t$narray[2]\t$rpm_mapped\t$narray[3]\t$narray[4]\t$exp\n";
				}
			}
		}
	}
}
chdir "../";
