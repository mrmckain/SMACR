#!/bin/bash

cd `pwd`

perl batch_filter_with_gff3.pl test.gff ../dmel-all-miRNA-r6.01.fasta ../dmel-all-ncRNA-r6.01.fasta ../dmel-all-transposon-r6.01.fasta
