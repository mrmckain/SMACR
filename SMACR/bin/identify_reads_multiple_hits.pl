#!/usr/bin/perl -w
use strict;

my %seqs;
my %badlist;
my %counts;
my %data;
chdir "Feature_Annotated";
my @files = <FullSet_Mismatch*_Results_*>;

for my $file (@files){
	open my $tfile, "<", $file;
	$file =~ /FullSet_Mismatch(\d+)_Results_(\w+).txt/;
        my $feature = $2;
        my $mismatch=$1;
	while(<$tfile>){
		chomp;
		
		my @tarray = split /\s+/;
		unless(%seqs){
			$seqs{$mismatch}{$feature}{$tarray[4]}=1;
		}
		for my $fea2 (keys %{$seqs{$mismatch}}){
			if(exists $seqs{$mismatch}{$fea2}{$tarray[4]}){
				$seqs{$mismatch}{$feature}{$tarray[4]}=1;
				unless($feature =~ /miRNA/ && $fea2 =~ /miRNA/){
					if($feature ne $fea2){
						if ($feature eq "Transposons" || $fea2 eq "Transposons" && $feature eq "Leftovers" || $fea2 eq "Leftovers"){
							$seqs{$mismatch}{$feature}{$tarray[4]}=1;
						}
						else{
							push(@{$badlist{$mismatch}{$tarray[4]}}, $feature);
		
							$seqs{$mismatch}{$feature}{$tarray[4]}=1;
						}
					}
		

				}
			}
	
					my $revcomp = reverse($tarray[4]);
					$revcomp =~ tr/ATCGatcg/TAGCtagc/;
			if(exists $seqs{$mismatch}{$fea2}{$revcomp}){
				$seqs{$mismatch}{$feature}{$tarray[4]}=1;
				unless($feature =~ /miRNA/ && $fea2 =~ /miRNA/){
                                        if($feature ne $fea2){
						if ($feature eq "Transposons" || $fea2 eq "Transposons" && $feature eq "Leftovers" || $fea2 eq "Leftovers"){
                                                        $seqs{$mismatch}{$feature}{$tarray[4]}=1;
                                                }
                                                else{
                                                        push(@{$badlist{$mismatch}{$tarray[4]}}, $feature);

                                                        $seqs{$mismatch}{$feature}{$tarray[4]}=1;
                                                }
					}
				}
			}
			$seqs{$mismatch}{$feature}{$tarray[4]}=1;
			
		}
		$seqs{$mismatch}{$feature}{$tarray[4]}=1;
	}
}
#open my $badout, ">", "Results_Found_Multiple_Places.txt";
# print $badout "Chromosome\tStart\tSize\tDirection\tSeq\tHits\tReadsPerMillionMapped\tReadsPerMillion\tRPMPerHits\tExperiment\n";
#open my $otherout, ">", "Multi_Hit_Features.txt";
#open my $countfile, ">", "PerExperiment_Counts.txt";
for my $file (@files){
        open my $tfile, "<", $file;
	open my $tout, ">>", $file . "_stripped.txt";
       	print $tout "Chromosome\tStart\tSize\tDirection\tSeq\tHits\tReadsPerMillionMapped\tReadsPerMillion\tRPMPerHits\tExperiment\n";
	$file =~ /FullSet_Mismatch(\d+)_Results_(\w+).txt/;
        my $feature = $2;
        my $mismatch=$1;
	open my $badout, ">>", "Results_Found_Multiple_Places_$mismatch.txt";
	if(-s "Results_Found_Multiple_Places_$mismatch.txt" == 0){
		print $badout "Chromosome\tStart\tSize\tDirection\tSeq\tHits\tReadsPerMillionMapped\tReadsPerMillion\tRPMPerHits\tExperiment\n";
	}
	open my $otherout, ">>", "Multi_Hit_Features_$mismatch.txt";
	while(<$tfile>){
                chomp;
		my @tarray = split/\s+/;
		
		my $revcomp = reverse($tarray[4]);
                $revcomp =~ tr/ATCGatcg/TAGCtagc/;
		
		if (exists $badlist{$mismatch}{$tarray[4]} || exists $badlist{$mismatch}{$revcomp}){
			print $badout "$_\n";
			if(exists $badlist{$mismatch}{$tarray[4]}){
				my $temp = join(",", @{$badlist{$mismatch}{$tarray[4]}});
				print $otherout "$tarray[4]\t$file\t$temp\n";
			}
			else{
				 my $temp = join(",", @{$badlist{$mismatch}{$revcomp}});
				print $otherout "$tarray[4]\t$file\t$temp\n";
			}
		}
		else{
			print $tout "$_\n";
			my @tarray = split /\s+/;
			$file =~ /FullSet_Mismatch\d_Results_(\w+).txt/;
                	my $feature = $1;
		#	$counts{$mismatch}{$feature}{$tarray[9]}++;
			$data{$mismatch}{$feature}{$tarray[9]}{$tarray[4]}{"RPMM"}=$tarray[6];
        		$data{$mismatch}{$feature}{$tarray[9]}{$tarray[4]}{"RPM"}=$tarray[7];
		}
	}
} 			
my %cleaned;

for my $mismatch (keys %data){
for my $feature (keys %{$data{$mismatch}}){
        for my $exp (sort {$a cmp $b} keys %{$data{$mismatch}{$feature}}){
                for my $seq (keys %{$data{$mismatch}{$feature}{$exp}}){
                        $cleaned{$mismatch}{$feature}{$exp}{"RPMM"}+=$data{$mismatch}{$feature}{$exp}{$seq}{"RPMM"};
                        $cleaned{$mismatch}{$feature}{$exp}{"RPM"}+=$data{$mismatch}{$feature}{$exp}{$seq}{"RPM"};
                        $counts{$mismatch}{$feature}{$exp}++;
                }
        }
}
}
for my $mismatch (keys %counts){
	open my $countfile, ">", "PerExperiment_Counts_Mismatch$mismatch.txt";
	#print $countfile "Feature\tExperiment\tCounts\n";
	print $countfile "Feature\tExperiment\tReadsPerMillionMapped\tReadsPerMillion\tTotalCount\n";
for my $feature (keys %{$cleaned{$mismatch}}){
        for my $exp (sort {$a cmp $b} keys %{$cleaned{$mismatch}{$feature}}){
                        print $countfile "$feature\t$exp\t$cleaned{$mismatch}{$feature}{$exp}{RPMM}\t$cleaned{$mismatch}{$feature}{$exp}{RPM}\t$counts{$mismatch}{$feature}{$exp}\n";
        }
}
}

chdir "../";
