#!/bin/bash

#Propreccessing the paired-end reads if it is necessary
python cut_reads.py reads_1.trimmed.fastq reads_2.trimmed.fastq reads_1.trimmed_cut.fastq reads_2.trimmed_cut.fastq

#Inverting the head of the fasta file in contig before mapping
awk 'BEGIN{id=1}{if($0~/>/){printf(">%d\n",id);id++}else{print $0}}' input_contigs.fa > contigs.fa

#Mapping the paired-end reads to the contigs using Bowtie2
bowtie2-build contigs.fa contigs
bowtie2 -x contigs -1 reads_1.trimmed_cut.fastq -2 reads_2.trimmed_cut.fastq -S contigs_short.sam

#iIverting the sam file to the bam file using samtools
samtools view -bS contigs_short.sam > contigs_short.bam
samtools sort contigs_short.bam contigs_short.sort
samtools index contigs_short.sort.bam

#Using MEC to detect and correct the misassemblies
python mec.py -bam contigs_short.sort.bam -i contigs.fa -o contigs-corr.fa