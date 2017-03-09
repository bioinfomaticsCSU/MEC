#!/usr/bin/env python

import os,sys
import commands
import pysam
import re
import parse_sam as ps
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

def delblankline(infile, outfile):
    """ Delete blanklines of infile """
    
    infp = open(infile, "r")
    outfp = open(outfile, "w")
    for li in infp:
        if li.split():
            outfp.writelines(li)
    infp.close()
    outfp.close()


def sub_header(contigFile):
    
    """
    Function to simplify the header of the input contig file.

    Arguments:
    contigPath: the filename or the path  of contigfile
    outcontigFile: the simplified contig file

    """
    contigfile = open(contigFile)
    outcontigFile = contigFile + ".s.fa"
    out = open(outcontigFile,"w")
    num = 0
    flag = False
    for line in contigfile:
        if not line: break
        if re.findall(r'^>',line):
            out.write(">"+str(num)+"\n")
            num = num+1
        else:
            out.write(line)
    contigfile.close()
    out.close()

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
            seq[contigid] += line
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
    GCRate = 0.0
    for cid in contigid:
        GCRate += GCProportion(seq[cid])
    GCavgRate = GCRate/len(seq)
    for cid in contigid:
        assert len(dic_interval[cid])>=1
        if len(dic_interval[cid])==1:
            if NProportion(seq[cid][dic_interval[cid][0][0]:dic_interval[cid][0][1]+1])>0.7: continue
            out.write(">"+cid+"\n"+seq[cid][dic_interval[cid][0][0]:dic_interval[cid][0][1]+1]+"\n")
        else:
            interval = dic_interval[cid]
            interval_len = len(interval)
            length = len(seq[cid])
            intervals = []
            p0 = interval[0][0]
            p1 = interval[0][1]
            for i in xrange(interval_len-1):
                start = max(interval[i][1]-100,0)
                end = min(interval[i+1][0]+100,length-1)
                if GCProportion(seq[cid][start:end]) < GCavgRate-0.1 or GCProportion(seq[cid][start:end]) > GCavgRate+0.1:
                    pos = i+1
                    p1 = interval[pos][1]
                    if pos == interval_len-1:
                        intervals.append((p0,p1))
                else:
                    intervals.append((p0,p1))
                    pos= i+1
                    p0 = interval[pos][0]
                    p1 = interval[pos][1]
                    if pos == interval_len-1:
                        intervals.append((p0,p1))
            cnt = 0
            for item in intervals:  
                if NProportion(seq[cid][item[0]:item[1]+1])>0.7: continue
                out.write(">"+cid+"_"+str(cnt)+"\n"+seq[cid][item[0]:item[1]+1]+"\n")
                cnt = cnt+1
                split_num += 1
            split_num -= 1
    out.close()
    return split_num


def cur_file_dir():
    path = sys.path[0]
    if os.path.isdir(path):
        return path
    elif os.path.isfile(path):
        return os.path.dirname(path)
                
def main():
    # read command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Tool to identify and correct misassemblies in de novo assemblies')
    parser.add_argument('-i', metavar='infasta', type=str, help='input assembly fasta file')
    parser.add_argument('-o', metavar='outfasta', type=str, help='output the correct file name')
    #parser.add_argument('-bam', metavar='bam', type=str, default=None,help='index bam file for alignment')
    parser.add_argument('-m1', metavar='mate1', type=str,help='one read pair for pair end reads')
    parser.add_argument('-m2', metavar='mate2', type=str,help='another read pair for pair end reads')
    parser.add_argument('-dir', metavar='direction', type=str, default="fr",help='the alignment direction')
    parser.add_argument('-c', metavar ='mincovrate', type=float, default=7.0, help='It is used to be divided by the average of the fragment coverage.')
    parser.add_argument('-r', metavar ='minreadrate', type=float, default=2.0 , help='the minimum coverage rate of mapped paired-end reads concordantly to the contigs to merge the interval')
    parser.add_argument('-q', metavar ='minmapq', type=int,default=40, help='Minimum MapQ value, above which a read pair is included in calculating population statistics. Default value: 40')
    parser.add_argument('-f', metavar ='minfraction', type=float, default=0.7, help='Minimum fraction of read pairs with correct orientation to call support for the assembly. Default value: 0.7')
    parser.add_argument('-n', metavar ='minnum', type=int, default=2, help='Minimum fraction of read pairs with correct orientation to call support for the assembly. Default value: 0.7')
    parser.add_argument('-sspace', metavar ='sspace', type=bool, default=0, help='whether or not choose the sspace scaffolding tool')
    args = parser.parse_args()


    sub_header(args.i)
    run = commands.getstatusoutput
    cmd = "mv " + args.i+".s.fa" + " " + args.i
    s, o = run(cmd)

    outfile = args.i + ".d.fa"
    delblankline(args.i,outfile)
    run = commands.getstatusoutput
    cmd = "mv " + outfile + " " + args.i
    s, o = run(cmd)
    
    third_party_path=cur_file_dir().replace("src","third_party")
    samfile = ps.generate_alignment(args.i,args.m1,args.m2,args.dir,third_party_path)
    bamfile = ps.sam_to_bam(samfile,third_party_path)+".bam"
   
    sam = pysam.Samfile(bamfile, "rb" )
    contigid = sam.references
    interval = pb.get_all_interval(sam,args.q,args.c,args.r,args.f,args.n,args.sspace)

    '''
    contigid = []
    if args.bam != None:
        bamfile = args.bam
        samfile = pysam.Samfile(bamfile, "rb" )
        contigid = samfile.references
        interval = pb.get_all_interval(samfile,args.q,args.c,args.r,args.f,args.n,args.sspace)
    '''
    
    
    print "start parse contig file"
    seq = parse_contig(args.i)

    print "output all the intervals for each contig"
    interval_file = "./intervals.txt"
    f = open(interval_file,"w")
    for cid in contigid:
        f.write(">"+cid+"\n")
        for item in interval[cid]:
            f.write(str(item[0])+" "+str(item[1])+"\n")
      
    print "output contig file"
    split_num = output_contig(seq,interval,contigid,args.o)
    outfile = args.o + ".d.fa"
    delblankline(args.o,outfile)
    run = commands.getstatusoutput
    cmd = "mv " + outfile + " " + args.o
    s, o = run(cmd)
    
    print "split_num",split_num
    print "\n"

    
if __name__ == "__main__":
    main()

