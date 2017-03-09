#!/usr/bin/env python

import re
import commands
import os

def clean_index(indexname):
    """
        Function to clean all of the index file with bowtie2.

        Arguments:
        filename: the index file path or name

     """
    for f in [indexname+".1.bt2",indexname+".2.bt2",indexname+".3.bt2",indexname+".4.bt2",indexname+".rev.1.bt2",indexname+".rev.2.bt2"]:
        if os.path.exists(f):
            os.remove(f)
    return

def clean_sam(samname):
    """
        Function to clean all of the sam file with bowtie2.

        Arguments:
        filename: the sam file path or name

     """
    if os.path.exists(samname):
        os.remove(samname)
    return

def generate_alignment(contigs,read1Path,read2Path,dirction,bowtie2):
    
    """
        Function to run commands to use
    """
    run = commands.getstatusoutput
    indexname = contigs.split(".")[0]
    sam = contigs.split(".")[0]+".sam"

    cmd = bowtie2+"/bowtie2-2.1.0/bowtie2-build " + contigs + " " + indexname
    print "==="
    print cmd
    s, o = run(cmd)
    if s != 0:
        print("bowtie2-build faild. exiting....")
        print(o)
        exit()
    print("build index complete...")

    readsformat = ""
    if re.findall(r'.fa$',read1Path) or re.findall(r'.fasta$',read1Path):
        readsformat = "-f"
    cmd = bowtie2+"/bowtie2-2.1.0/bowtie2 " + readsformat + " -x " + indexname + " -1 " + read1Path + " -2 " + read2Path + " -S " + sam + " --" + dirction
    print cmd
    s, bt2output = run(cmd)
    print(bt2output)
    if s != 0:
        print("bowtie2 faild. exiting....")
        print(bt2output)
        exit()
    print("mapping complete...")
    print("====")
    clean_index(indexname)
    return sam
     
def sam_to_bam(samfile,samtool):

    """
        Function to convert sam file to bam file.
    """
    run = commands.getstatusoutput
    bam = samfile + ".bam"
    cmd = samtool+"/samtools-0.1.18/samtools view -bS " + samfile + " >" + bam
    print cmd
    s, o = run(cmd)
    if s != 0:
        print("samtools view faild. exiting....")
        print(o)
        exit()
    bam_sort = bam +".sort"
    cmd = samtool+"/samtools-0.1.18/samtools sort " + bam + " "+ bam_sort
    print cmd
    s, o = run(cmd)
    if s != 0:
        print("samtools sort faild. exiting....")
        print(o)
        exit()
    cmd =samtool+"/samtools-0.1.18/samtools index " + bam_sort+".bam"
    print cmd
    s, o = run(cmd)
    if s != 0:
        print("samtools index faild. exiting....")
        print(o)
        exit()
    clean_sam(samfile)
    return bam_sort

def concordant_rate(bowtie2_msg):
    patt = re.compile("\(([0-9.]*)%\) aligned concordantly exactly 1 time")
    res = patt.search(bowtie2_msg).groups()
    x = float(res[0])
    patt = re.compile("([0-9.]*)% overall alignment rate")
    res = patt.search(bowtie2_msg).groups()
    y = float(res[0])
    return x/y
