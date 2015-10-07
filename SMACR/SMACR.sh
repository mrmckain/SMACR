#!/bin/bash

####################################################
#Sequence Mapping, Annotation, and Counting for sRNAS
#
# Pipeline for analyzing small RNA Illumina reads on the Drosophila genome

#USAGE 1=DATA DIRECTORY

TIME=`date "+%m_%d_%y"`
mkdir SMACR_Run_$TIME
cd SMACR_Run_$TIME
mkdir Cleaned_Reads
cd Cleaned_Reads

for f in $1/*fastq.gz
do
	echo $f
	EXPNAME=$(echo $f | perl -ne 'print $1 if /\/\/(.*).fastq.gz/'| cut -d "_" -f 1)
	#EXPNAME=$(echo $f | cut -d "_" -f 1)
	java -jar /data/steinigerlab/sRNAcleanreads2/SMACR/bin/Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 $f $EXPNAME\_trimmed.fq ILLUMINACLIP:/data/steinigerlab/sRNAcleanreads2/SMACR/bin/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:1:20:7 SLIDINGWINDOW:3:20 MINLEN:19
	perl /data/steinigerlab/sRNAcleanreads2/SMACR/bin/max_read_filter.pl $EXPNAME\_trimmed.fq
done
cd ../

perl /data/steinigerlab/sRNAcleanreads2/SMACR/bin/srna_processing.pl

perl /data/steinigerlab/sRNAcleanreads2/SMACR/bin/mapping_read_stats.pl

perl /data/steinigerlab/sRNAcleanreads2/SMACR/bin/reformat_srna_output.pl

perl /data/steinigerlab/sRNAcleanreads2/SMACR/bin/filter_results_byexperiment.pl

perl /data/steinigerlab/sRNAcleanreads2/SMACR/bin/batch_filter_with_gff3.pl /data/steinigerlab/sRNAcleanreads2/SMACR/genome_anno/dmel-all-r6.01.gff /data/steinigerlab/sRNAcleanreads2/SMACR/genome_anno/dmel-all-miRNA-r6.01.fasta /data/steinigerlab/sRNAcleanreads2/SMACR/genome_anno/dmel-all-ncRNA-r6.01.fasta /data/steinigerlab/sRNAcleanreads2/SMACR/genome_anno/dmel-all-transposon-r6.01.fasta

perl /data/steinigerlab/sRNAcleanreads2/SMACR/bin/identify_reads_multiple_hits.pl

perl /data/steinigerlab/sRNAcleanreads2/SMACR/bin/raw_rpmm_srna_size.pl

perl /data/steinigerlab/sRNAcleanreads2/SMACR/bin/raw_rpmm_srna_bases.pl 
