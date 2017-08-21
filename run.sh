#!/bin/bash

#In the following, there is an example of command lines with MEC to correct misassemblies.

#Propreccessing the paired-end reads
python cut_reads.py reads_1.trimmed.fastq reads_2.trimmed.fastq reads_1.trimmed_cut.fastq reads_2.trimmed_cut.fastq

#Inverting the head of the fasta file in contig before mapping
awk 'BEGIN{id=1}{if($0~/>/){printf(">%d\n",id);id++}else{print $0}}' input_contigs.fa > contigs.fa

#Mapping the paired-end reads to the contigs using Bowtie2
bowtie2-build contigs.fa contigs
bowtie2 -x contigs -1 reads_1.trimmed_cut.fastq -2 reads_2.trimmed_cut.fastq -S contigs_short.sam

#Inverting the sam file to the bam file using samtools
samtools view -bS contigs_short.sam > contigs_short.bam
samtools sort contigs_short.bam contigs_short.sort
samtools index contigs_short.sort.bam

#Deleting the unneccessary files
rm -rf *.bt2
rm contigs_short.sam
rm contigs_short.bam

#Using MEC to detect and correct the misassemblies
python mec.py -bam contigs_short.sort.bam -a 0.4 -b 0.5 -m 600 -s 100 -i contigs.fa -o contigs-corr.fa
