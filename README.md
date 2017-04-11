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

	If the users don't install python and pysam module in your PATH, you must first install it. 
	Or you can use the executable python in the directory /third_party/python2.7/bin/python, which is packaged with the source code, and put it in your PATH.

3) Running
```
	Please go to the directory "src".
	Run command line:  
	python mec.py -m1 pair1.fasta -m2 pair2.fasta -i assembly.fasta -o correct_assembly.fasta [options] 
	[options]
	-i <input_assembly.fasta>
		Mandatory parameter. The input assembly fasta file.
	-o <output_correct_assembly.fasta>
		Mandatory parameter. The output corrected fasta file.
	-m1 <mate1>
		Mandatory parameter. The mate1 of paired-end reads mapped to the contigs.
	-m2 <mate2>
		Mandatory parameter. The mate2 paired-end reads mapped to the contigs.
	-dir
		Optional parameter. The direction of the pair end reads. Default value: fr.
	-c <minimum fragment coverage rate>
		Optional parameter. It is used to be divided by the average of the fragment coverage. Default value: 7.0.
	-r <minimum read coverage rate>
		Optional parameter. It is used to be decide the minmum reads coverage to remove false misassemblies. Default value: 2.0.
	-q <minimum mapping quality>
		Optional parameter. Only the pair end reads whose mapping quality larger than minimum mapping quality to compute the fragment coverage. Default value: 40.
	-f <minimum fraction>
		Optional parameter. Minimum fraction of read pairs whose insert size distrbution is accordant and discordant to validate the correct assembly. Default value: 0.7.
	-n <minimum number>
		Optional parameter. Minimum number of read pairs which map discordantly on the candidate misassembly. Default value: 2.
	-sspace
		Optional parameter. Whether use sspace to scaffold. Default value: 0.
```	
4) Output:

	The final output file including the corrected fasta file and the correct interval for each contig ("intervals.txt").

	
	
	
	
	
