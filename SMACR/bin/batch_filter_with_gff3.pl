#!/usr/bin/env perl -w
use strict;
#USAGE:perl batch_filter_with_gff3.pl full_gff3_file miRNA_sequence_file ncRNA_sequence_file transposon_sequence_file

`mkdir Feature_Annotated`;
chdir "Feature_Annotated";
my %targets;
my $miRNA_file = $ARGV[1];
my $ncRNA_file = $ARGV[2];
my $transposons_file = $ARGV[3];

my %feature_ids;

open my $tfile, "<", $miRNA_file;
while(<$tfile>){
	chomp;
	if(/^>/){
		if(/pre_miRNA/){
			my $tid = substr($_, 1, index($_, " ")-1);
                	$feature_ids{"pre_miRNA"}{$tid}=1;
		}
		else{
			my $tid = substr($_, 1, index($_, " ")-1);
			$feature_ids{"miRNA"}{$tid}=1;
		}
	}
}

open $tfile, "<", $ncRNA_file;
while(<$tfile>){
        chomp;
        if(/^>/){
                my $tid = substr($_, 1, index($_, " ")-1);
                $feature_ids{"ncRNA"}{$tid}=1;
        }
}

open $tfile, "<", $transposons_file;
while(<$tfile>){
        chomp;
        if(/^>/){
                my $tid = substr($_, 1, index($_, " ")-1);
                $feature_ids{"Transposons"}{$tid}=1;
        }
}


open my $gff, "<", $ARGV[0]; #gff file
while(<$gff>){
	next if (/^\#/);
	chomp;
	if(/ID=(FB\w+)/){
		my $tid=$1;
	if($tid){
		for my $type (keys %feature_ids){
			if(exists $feature_ids{$type}{$tid}){
				
				my @tarray = split /\s+/;

				$targets{$tarray[0]}{$type}{$tarray[3]}{$tarray[4]}=1;
			}
		}
	}}
}

#harcoding the hairpin stuff on X and 2L

$targets{"X"}{"Hairpin"}{4921600}{4927700}=1;
$targets{"2L"}{"Hairpin"}{9788840}{9790180}=1;

my @files = <../Experiment_Filtered_Results/*Results_ExpFiltered.txt>;

for my $file (@files){

	my %found_sites;	

	my ($region, $match);
        $file =~ /Results\/(.+)\_(Mismatch\d).+/;
        $region=$1;
        $match=$2;
#adde the hairpin seqs, too
#	open my $miRNA_split_file, ">", $region . "_" . $match . "_Results_" . "miRNA" . ".txt"; #entry line is type to split out
	open my $no_split, ">", $region . "_" . $match . "_Results_Leftovers.txt";
=item
	open my $miRNA_trans_split, ">", $region . "_" . $match . "_Results_TransSplit_" . "miRNA" . ".txt";
	open my $ncRNA_split_file, ">", $region . "_" . $match . "_Results_" . "ncRNA" . ".txt";	
	open my $Transposons_split_file, ">", $region . "_" . $match . "_Results_" . "Transposons" . ".txt";
	open my $ncRNA_trans_split, ">", $region . "_" . $match . "_Results_TransSplit_" . "ncRNA" . ".txt";
	open my $Transposons_trans_split, ">", $region . "_" . $match . "_Results_TransSplit_" . "Transposons" . ".txt";
	open my $Hairpin_trans_split, ">", $region . "_" . $match . "_Results_TransSplit_" . "Hairpin" . ".txt";
	open my $Hairpin_split_file, ">", $region . "_" . $match . "_Results_" . "Hairpin" . ".txt";  

=cut

	open my $tfile, "<", $file;
	
	my $header = readline($tfile);
=item
	print $miRNA_split_file "$header";
	print $no_split "$header";
	print $miRNA_trans_split "$header";
	print $ncRNA_split_file "$header";
        print $ncRNA_trans_split "$header";
	print $Transposons_split_file "$header";
        print $Transposons_trans_split "$header";
	print $Hairpin_split_file "$header";
	print $Hairpin_trans_split "$header";
=cut
	my $found_it;
 	#SRNA:	
	while(<$tfile>){
		chomp;
		$found_it=();
		my @tarray = split /\s+/;
		if ($tarray[0] =~ /arm/){
			$tarray[0] =~ /(arm)_(.+)/;
			$tarray[0] = $2;
		}
		my $temp_type;
		my $temp_start;
		my $temp_end;
		for my $type (keys %{$targets{$tarray[0]}}){
		for my $start (sort {$a<=>$b} keys %{$targets{$tarray[0]}{$type}}){
		for my $end (keys %{$targets{$tarray[0]}{$type}{$start}}){
			if($temp_start){
				if ($start == $temp_start){
					if($end == $temp_end){
						next;
					}
				}
			}
			if ($tarray[1] >= $start){
				if($tarray[1] <= $end){
					if(($tarray[1]+$tarray[2]-1) <= $end){
				#		print "$type\t$start\t$end\t$_\n";
						open my $tfilename, ">>",  $region . "_" . $match . "_Results_" . $type . ".txt";
						my $filesize = -s $region . "_" . $match . "_Results_" . $type . ".txt";
						if($filesize == 0){
							print $tfilename "$header";
						}
						print $tfilename "$_\n";
						$found_it=1;
						$temp_type = $type;
						$temp_start = $start;
						$temp_end = $end;
						#goto SRNA;
					}
					else{
						#print "$type\t$start\t$end\t$_\n";
						open my $tfilename, ">>",  $region . "_" . $match . "_Results_TransSplit_" . $type . ".txt";
						my $filesize = -s $region . "_" . $match . "_Results_TransSplit_" . $type . ".txt";
                                                if($filesize == 0){
                                                        print $tfilename "$header";
                                                }
						print $tfilename "$_\n";
						$found_it=1;
						$temp_type = $type;
                                                $temp_start = $start;
                                                $temp_end = $end;
						#goto SRNA;
					}
				}
			}
			elsif($tarray[1] < $start && ($tarray[1]+$tarray[2]-1) > $start && ($tarray[1]+$tarray[2]-1) <= $end){
				#print "$type\t$start\t$end\t$_\n";
				open my $tfilename, ">>",  $region . "_" . $match . "_Results_TransSplit_" . $type . ".txt";
                       		my $filesize = -s $region . "_" . $match . "_Results_TransSplit_" . $type . ".txt";
                                if($filesize == 0){
                                	print $tfilename "$header";
                                }
                               	print $tfilename "$_\n";
                              	$found_it=1;
				$temp_type = $type;
                                                $temp_start = $start;
                                                $temp_end = $end;
                                 #goto SRNA;
			}
		}
		}}
		unless ($found_it){
			print $no_split "$_\n";
			$temp_type = ();
                                                $temp_start = ();
                                                $temp_end = ();
		}
	}
}		 		
