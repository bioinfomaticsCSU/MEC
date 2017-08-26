License
=========

Copyright (C) 2014 Jianxin Wang(jxwang@mail.csu.edu.cn), Binbin Wu(binbinwu@csu.edu.cn)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.

Jianxin Wang(jxwang@mail.csu.edu.cn), Binbin Wu(binbinwu@csu.edu.cn)
School of Information Science and Engineering
Central South University
ChangSha
CHINA, 410083


MEC
=================
1) Introduction

	MEC is a misassembly identify and correction tool.
	The input includes the paired-end read and original contigs. 

2) Installing

	First, please build and install Python and pysam module. 
	The version of python we used is Python 2.7.11, which can be downloaded from https://www.python.org/downloads/release/python-2711/.
	The version of Cython  we used is Cython-0.25.2, which can be downloaded from https://pypi.python.org/pypi/Cython/0.25.2.
	The version of pysam module we used is pysam-0.8.4, which can be downloaded from https://pypi.python.org/pypi/pysam/0.8.4. 
	
	In addition, please build and install bowtie2 and samtools.
	The version of Bowtie2 we used is bowtie2-2.2.8, which can be downloaded from https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2.
	The version of Samtools we used is samtools-1.2, which can be downloaded from https://sourceforge.net/projects/samtools/files/samtools/.
	
3) Running
```
	Please go to the directory "src".
	Run command line:  
	python mec.py -bam pairs.sort.bam -i assembly.fasta -o correct_assembly.fasta [options] 
	[options]
	-i <input_assembly.fasta>
		Mandatory parameter. The input assembly fasta file.
	-o <output_correct_assembly.fasta>
		Mandatory parameter. The output corrected fasta file.
	-bam <the index bam file>
		Mandatory parameter. The index bam file for alignment. 
	-q <minimum mapping quality>
		Optional parameter. The minimum value of mapping quality. Default value: 40.
	-m <mu>
		Optional parameter. The mean of the paired-end reads' insert size. Default value: 0.
	-s <sigma>
		Optional parameter. The variance of the paired-end reads' insert size. Default value: 0.
	-a <alpha>
		Optional parameter. The percentage of the average of the fragment coverage. Default value: 0.4.
	-b <beta>
		Optional parameter. One cutoff for removing false misassemblies. Default value: 0.5.
	-g <gamma>
		Optional parameter. One parameter for determining high or low GC content. Default value: 0.2.
```	

4) Running
```
	Please go to the directory "src".
	Run command line:  
	python mec.py -bam pairs.sort.bam -i assembly.fasta -o correct_assembly.fasta [options] 
	[options]
	-i <input_assembly.fasta>
		Mandatory parameter. The input assembly fasta file.
	-o <output_correct_assembly.fasta>
		Mandatory parameter. The output corrected fasta file.
	-bam <the index bam file>
		Mandatory parameter. The index bam file for alignment. 
	-q <minimum mapping quality>
		Optional parameter. The minimum value of mapping quality. Default value: 40.
	-m <mu>
		Optional parameter. The mean of the paired-end reads' insert size. Default value: 0.
	-s <sigma>
		Optional parameter. The variance of the paired-end reads' insert size. Default value: 0.
	-a <alpha>
		Optional parameter. The percentage of the average of the fragment coverage. Default value: 0.4.
	-b <beta>
		Optional parameter. One cutoff for removing false misassemblies. Default value: 0.5.
	-g <gamma>
		Optional parameter. One parameter for determining high or low GC content. Default value: 0.2.
```	

5) Output:

	The final output file including the corrected fasta file and the correct interval for each contig ("intervals.txt").
