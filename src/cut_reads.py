#! /usr/bin/python

import os,sys

def cut_reads01(file01,file02):

    n =0
    inf = open(file01)
    outf = open(file02,"w")
    for line in inf:
        l = line.replace("\n","")
        n+=1
        if n%4==1:
            outf.write(l+"\n")
        elif n%4==2:
            if len(l)>=100:
                outf.write(l[0:100]+"\n")
            else:
                outf.write(l+"\n")
        elif n%4==3:
            outf.write(l+"\n")
        else:
            if len(l)>=100:
                outf.write(l[0:100]+"\n")
            else:
                outf.write(l+"\n")
            
    inf.close()
    outf.close()

def cut_reads02(file01,file02):
    n =0
    inf = open(file01)
    outf = open(file02,"w")
    for line in inf:
        l = line.replace("\n","")
        n+=1
        if n%4==1:
            outf.write(l+"\n")
        elif n%4==2:
            if len(l)>=100:
                #outf.write(l[len(l)-100:]+"\n")
                outf.write(l[0:100]+"\n")
            else:
                outf.write(l+"\n")
        elif n%4==3:
            outf.write(l+"\n")
        else:
            if len(l)>=100:
                #outf.write(l[len(l)-100:]+"\n")
                outf.write(l[0:100]+"\n")
            else:
                outf.write(l+"\n")
    inf.close()
    outf.close()

def main():
    # read command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Tool to cut the paired-end reads')
    parser.add_argument('r1', metavar='read1', type=str, help='the input left reads of the paired-end reads')
    parser.add_argument('r2', metavar='read2', type=str, help='the input right reads of the paired-end reads')
    parser.add_argument('o1', metavar='out1', type=str, help='the output of left reads of the paired-end reads')
    parser.add_argument('o2', metavar='out2', type=str, help='the output of right reads of the paired-end reads')
    args = parser.parse_args()

    file01 = args.r1
    out01 = args.o1
    file02 = args.r2
    out02 = args.o2
    cut_reads01(file01,out01)
    cut_reads02(file02,out02)
    
    
if __name__ == "__main__":
    main()
