#!/usr/bin/perl -w
`mkdir sRNA_Mapping`;
chdir "sRNA_Mapping";
my @files = <../Cleaned_Reads/*reads.cleaned.fq>;

for my $rfile (@files){

	my %srnas;

	open my $reads, "<", $rfile;
	while(<$reads>){

    		my $seq = readline($reads);
        	chomp($seq);
		my $rseq = reverse($seq);
		$rseq =~ tr/ATCGatcg/TAGCtagc/;
		if(exists $srnas{$rseq}){
			$srnas{$rseq}++;
		}
		else{
			$srnas{$rseq}++;
		}
        	readline($reads);
        	readline($reads);
	}


	my $count = keys %srnas;
	$rfile =~ /\/(\w+)_reads.cleaned.fq/;
	my $experiment = $1;
	my $outfilename = $experiment . ".srnas.fsa";
	
	open my $out, ">", $outfilename;

	for my $ss (sort keys %srnas){
    		print $out ">$srnas{$ss}\n$ss\n";
	}

	for (my $mismatch = 0; $mismatch <= 1; $mismatch++){
		my $bowtieout = $experiment . "_Dmel6.0.1_mismatch$mismatch.sam";
		my $alf = $experiment . "_aligned_mismatch$mismatch" . ".fsa";
		my $unf = $experiment . "_unaligned_mismatch$mismatch" . ".fsa";
		`bowtie -a --al $alf --un $unf -v $mismatch -S -f /data/steinigerlab/sRNAcleanreads2/SMACR/genome_anno/Dmel6.0.1 $outfilename > $bowtieout`;
		
		my %hits;
		open my $out, ">", $experiment . "_sRNA_Dmel6.0.1_Mismatch$mismatch" . "_Results.txt";
		open my $file, "<", $bowtieout; #sam file

		my %true_hit_count;
		while(<$file>){
        		chomp;
        		if(/\@SQ/ || /\@HD/ || /\@PG/){
                		next;
        		}
       			my @tarray = split /\t/;
        
        		if($tarray[1] == 4){
                		next;
        		}
        		my $size = length($tarray[9]);
        		my $dir="+";
        		if($tarray[1] ==16){
                		$dir ="-";
        		}
        		if($dir eq "-"){
                		$tarray[9] =~ tr/ATCGatcg/TAGCtagc/;
                		$tarray[9] = reverse($tarray[9]);
        		}
        		$true_hit_count{$tarray[9]}++;
        		$hits{$tarray[2]}{$tarray[3]}{$size}{$dir}{$tarray[9]}++;

		}

		open my $reads, "<", $outfilename; #reads file

		my %reads;
		my $count;
		my $total_count=0;
		while(<$reads>){
        		chomp;
        		if (/>/){
                		$count = substr($_,1);
                		$total_count += $count;
        		}
        
        		else{
                		$reads{$_}=$count;
        		}
		}
		my $total_readsmill = $total_count/(1000000);

		print $out "Chromosome\tStart\tSize\tDirection\tSeq\tHits\tReadsPerMillion\tRPMPerHits\tExperiment\n"; 

		for my $chr (sort keys %hits){
        		for my $start ( sort {$a <=> $b} keys %{$hits{$chr}}){
                		for my $siz (sort {$a <=> $b} keys %{$hits{$chr}{$start}}){
                        		for my $direction (keys %{$hits{$chr}{$start}{$siz}}){
                                		for my $seq (keys %{$hits{$chr}{$start}{$siz}{$direction}}){
                                        		my $rpm = $reads{$seq}/$total_readsmill;
                                        		my $rpmph = $rpm/$true_hit_count{$seq};
                                        		print $out "$chr\t$start\t$siz\t$direction\t$seq\t$true_hit_count{$seq}\t$rpm\t$rpmph\t$experiment\n";
                                		}		
                        		}
                		}
        		}
		}
	}
}


chdir "../";

