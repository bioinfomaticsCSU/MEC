#!/usr/bin/env python

import os,sys
import pysam
import re
import parse_bam as pb

def NProportion(line):
    c = line.count("n") + line.count("N")
    t = len(line)*1.0
    if t==0: return 0
    return c/t

def GCProportion(line):
    c = line.count("g") + line.count("c")+line.count("G") + line.count("C")
    t = len(line)*1.0
    if t==0: return 0
    return c/t

def parse_contig(contigPath):

    """
    Function to parse the input contig file.
    Return a dict of (contigid:sequence).

    Arguments:
    contigPath: the filename or the path  of contigfile.
    flag: the choice to return seq or contigid.

    """
    contigfile = open(contigPath)
    seq={}
    contigid=""
    for line in contigfile:
        if not line: break
        if re.findall(r'^>',line):
            contigid = line.split(" ")[0][1:-1]
            seq[contigid] = ""
        else:
            seq[contigid] += line.strip()
    contigfile.close()
    return seq

def output_contig(seq,dic_interval,contigid,outfile):
    
    """
    Function to split the contigs at positions identified as assembly errors and write a new fasta file containing all contigs.
    
    Arguments:
    outfile: name of the new fasta file (including filepath).
    interval: list of range that have no misassemblies errors.
    contigid: list of id of the all contig.
    seq: the dictionary of the all contigs, key=id of the contig, value=the seq of each contig.
    """  
    out = open(outfile, "w")
    split_num = 0
    
    for cid in contigid:
        if len(dic_interval[cid])<1:print cid
        assert len(dic_interval[cid])>=1
        if len(dic_interval[cid])==1:
            p0=dic_interval[cid][0][0]
            p1=dic_interval[cid][0][1]
            if NProportion(seq[cid][p0:p1+1])>0.7: continue
            out.write(">"+cid+"\n"+seq[cid][p0:p1+1]+"\n")
        else:
            intervals = dic_interval[cid]
            cnt = 0
            for item in intervals:  
                if NProportion(seq[cid][item[0]:item[1]+1])>0.7: continue
                out.write(">"+cid+"_"+str(cnt)+"\n"+seq[cid][item[0]:item[1]+1]+"\n")
                cnt = cnt+1
            split_num += len(intervals)-1
    out.close()
    return split_num
                
def main():
    # read command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Tool to identify and correct misassemblies in de novo assemblies')
    parser.add_argument('-i', metavar='infasta', type=str, help='input assembly fasta file')
    parser.add_argument('-o', metavar='outfasta', type=str, help='output the correct file name')
    parser.add_argument('-bam', metavar='bam', type=str, default=None,help='index bam file for alignment')
    parser.add_argument('-q', metavar ='minmapq', type=int,default=40, help='Minimum mapping quality value. Default value: 40')
    parser.add_argument('-a', metavar ='alpha', type=float, default=0.4, help='The percentage of the average of the fragment coverage. Default value: 0.4')
    parser.add_argument('-b', metavar ='beta', type=float, default=0.5, help='One cutoff for removing false misassemblies. Default value: 0.5')
    parser.add_argument('-g', metavar ='gamma', type=float, default=0.15, help='One parameter for determining high or low GC content. Default value: 0.15.')    
    args = parser.parse_args()

    print "start parse contig file"
    seq = parse_contig(args.i)
    
    if args.bam != None:
        bamfile = args.bam
        samfile = pysam.Samfile(bamfile, "rb" )
        contigid = samfile.references
        GCRate = 0.0
        for cid in contigid:
            GCRate += GCProportion(seq[cid])
        GCavgRate = GCRate/len(seq)
        interval = pb.get_all_interval(samfile,args.q,args.a,args.b,args.g,seq,GCavgRate)

    else:
        print "error"
        
    print "output all the intervals for each contig"
    interval_file01 = "./intervals.txt"
    f01 = open(interval_file01,"w")
    for cid in contigid:
        f01.write(">"+cid+"\n")
        for item in interval[cid]:
            f01.write(str(item[0])+" "+str(item[1])+"\n")
    f01.close()
      
    print "output contig file"
    split_num = output_contig(seq,interval,contigid,args.o)
    
    print "print the number of detected misassemblies: ",split_num
    print "\n"
    
if __name__ == "__main__":
    main()

