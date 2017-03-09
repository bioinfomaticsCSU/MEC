#!/usr/bin/env python

import pysam
import collections
import re
import numpy as np

def get_read_size_distribution(all_reads,minmapq):

    """
    Function to calculate global insert size distribution across the whole assembly
    Return a frequency table of insert sizes as a  dictionary with key = insert size, value = frequency

    """

    frq = collections.defaultdict(int) # dictionary of insert sizes
    found = {}
    for read in all_reads:
        # accept read based on mapq, contig alignemnt and insert size
        if (read.mapq > minmapq) and (read.rnext == read.tid):
            if read.qname in found and found[read.qname][0]==read.tid:
                mate = found[read.qname]
                isize = abs(max( mate[1]+mate[2]-read.pos,read.pos+read.rlen-mate[1]))
                frq[isize] += 1
            else:
                found[read.qname] = (read.tid,read.pos,read.rlen)
    return frq

def MAD(frq):

    """

    Function to calculate median and median absolute deviation (MAD) from a dictionary of frequencies.

    Return a tuple of (mean, MAD).

    Arguments:
    frq: a dictionary of frequencies: key = insert size, value = frequency
     
    """

    all_lengths = []
    for k, v in frq.iteritems():
        new_vals = [k] * int(v)
        all_lengths.extend(new_vals)
    all_lengths = np.array(sorted(all_lengths))
    mid = len(all_lengths)/2
    median = all_lengths[mid]
    residuals = sorted(abs(all_lengths - median)) # difference between val and median
    MAD = residuals[mid] # median of residuals
    isize_sd = int(1.4826 * MAD)
    return median, isize_sd

def get_all_interval(samfile,minq,minrate,mincov,fraction,minnum,sspace):
    
    """
        Function to get all of the interval containing no misjoin errors among all of contigs.

        Arguments:
        filename: the bam file path or name
     """
    # getting reads from bamfile
    references = samfile.references  #type:list, item type:str
    lengths = samfile.lengths

    all_reads = samfile.fetch()
    frq = get_read_size_distribution(all_reads,minq)
    s,d = MAD(frq)
    #print s,d
    #print minrate,mincov,fraction
    interval = {}
    interval_valid = []
    for k,v in zip(references,lengths):
        reads = samfile.fetch(k, 0, v)
        interval_valid = get_valid_interval(k,0,v,lengths,reads,s,d,minq,minrate,mincov,fraction,minnum,sspace)
        interval[k] = interval_valid
    return interval
    
def get_valid_interval(ref,start,end,lengths,reads,mu,sigma,minq,minrate,mincov,fraction,minnum,sspace):

    cid=ref
    length = end
    left_accordant = []
    left_discordant = []
    coverage = [0]*length
    read_coverage=[0]*length
    
    for read in reads:
        if read.cigarstring!=None:
            if read.rnext == read.tid:
                for i in xrange(read.pos,min(read.pos+read.query_alignment_length,length)):
                    read_coverage[i]+=1
                if read.is_proper_pair:
                    if read.mapq>=minq:
                        for i in xrange(read.pos+read.rlen,read.pnext):
                            coverage[i]+=1
                if max(0,mu-3*sigma) <= abs(read.isize) <= mu+3*sigma:
                    if not read.is_reverse and read.mate_is_reverse:
                        if read.pnext>read.pos:
                            left_accordant.append((read.pos,read.pnext,read.isize))
                        else:
                            left_discordant.append((read.pos,read.pnext,read.isize))
                    elif read.is_reverse and not read.mate_is_reverse:continue
                    elif read.is_reverse == read.mate_is_reverse:
                        if read.pnext>read.pos:
                            left_discordant.append((read.pos,read.pnext,read.isize))
			if read.pnext<read.pos:
			    left_discordant.append((read.pnext,read.pos,read.isize))
                    else:
                        print read
                        print "else",read.isize,read.is_reverse,read.mate_is_reverse
                else:
                    if not read.is_reverse and read.mate_is_reverse:
                        if read.pnext>read.pos:
                            left_discordant.append((read.pos,read.pnext,read.isize))
			if read.pnext<read.pos:
			    left_discordant.append((read.pnext,read.pos,read.isize))
                    elif read.is_reverse and not read.mate_is_reverse:
                        continue
                    elif read.mate_is_reverse == read.is_reverse:
                        if read.pnext>read.pos:
                            left_discordant.append((read.pos,read.pnext,read.isize))
			if read.pnext<read.pos:
			    left_discordant.append((read.pnext,read.pos,read.isize))
                    else:
                        print read
                        print read.isize,read.is_reverse,read.mate_is_reverse
            else:
                if min(read.pos,length-(mu-3*sigma))+min(read.pnext,lengths[read.rnext]-(mu-3*sigma))>mu+3*sigma or min(read.pos,length-(mu-3*sigma))+min(read.pnext,lengths[read.rnext]-(mu-3*sigma))<mu-3*sigma:
                    left_discordant.append((read.pos,-1,read.isize))

    threshold = sum(coverage)/length
    adjacent_pos = [i for i in xrange(len(coverage)) if coverage[i]*minrate > threshold]
    frg_interval = pos_to_interval(adjacent_pos)

    read_threshold = sum(read_coverage)/length

    if len(frg_interval)>=2:
        interval_len = len(frg_interval)
        flag=True
        interval=[]
        p0 = frg_interval[0][0]
        p1 = frg_interval[0][1]
        for k in xrange(p0-1,-1,-1):
            if read_coverage[k]*minrate>=read_threshold:
                p0-=1
            else:
                break
        for i in xrange(interval_len-1):
            for j in xrange(frg_interval[i][1]+1,frg_interval[i+1][0]+1):
                if read_coverage[j]*mincov>=read_threshold:
                    p1+=1
                    flag=True
                else:
                    flag=False
                    break
            if flag==True:
                pos = i+1
                p1 = frg_interval[pos][1]
                if pos == interval_len-1:
                    for k in xrange(p1+1,length):
                        if read_coverage[k]*minrate>=read_threshold:
                            p1+=1
                        else:
                            break
                    interval.append((p0,p1))
            else:
                interval.append((p0,p1))
                pos= i+1
                p0 = frg_interval[pos][0]
                p1 = frg_interval[pos][1]
                if pos == interval_len-1:
                    for k in xrange(p1+1,length):
                        if read_coverage[k]*minrate>=read_threshold:
                            p1+=1
                        else:
                            break
                    interval.append((p0,p1))  
    else:
        if sspace==0:
            interval=[(0,length)]
        else:
            read_adjacent_pos = [i for i in xrange(len(read_coverage)) if read_coverage[i]*minrate>read_threshold]
            read_interval = pos_to_interval(read_adjacent_pos)
            if len(read_interval)>=1:
                interval=[(read_interval[0][0],read_interval[-1][1])]
            else:
                interval=[(0,length)] 

    if len(interval)>=2:
        if length >= mu-3*sigma:
            interval_len = len(interval)
            interval_valid = []
            p0 = interval[0][0]
            p1 = interval[0][1]
            for i in xrange(interval_len-1):
                start = max(interval[i][1]-mu,0)
                end = min(interval[i+1][0]+mu,length-1)
                isize_discordant = [item for item in left_discordant if ( (item[0]>=start and item[0]<=interval[i+1][0])  or (item[1]>=interval[i][1] and item[1]<=end) )]
                isize_accordant = [item for item in left_accordant if ( (item[0]>=start and item[0]<=interval[i+1][0]) and (item[1]>=interval[i][1] and item[1]<=end) )]  #or (item[1]>=interval[i][1] and item[1]<=end)
                accordante_rate = len(isize_accordant)*1.0/max((len(isize_accordant)+len(isize_discordant)),0.01) 
                if accordante_rate >= fraction or len(isize_accordant)<=minnum:
                    pos = i+1
                    p1 = interval[pos][1]
                    if pos == interval_len-1:
                        interval_valid.append((p0,p1))
                else:
                    interval_valid.append((p0,p1))
                    pos= i+1
                    p0 = interval[pos][0]
                    p1 = interval[pos][1]
                    if pos == interval_len-1:
                        interval_valid.append((p0,p1))
        else:
            interval_valid=interval
    else:
        interval_valid=interval
        
    return interval_valid
                
def pos_to_interval(pos_list):
    
    """
        Function to invert pos list to the interval if pos is contiguous.
        return: the interval that does not contain misassembly error.

        Arguments:
        pos_list: the list of pos

     """
    interval = []
    if len(pos_list)>=2:
        start = pos_list[0]
        end = pos_list[0]
        for i in pos_list[1:]:
            if i-end == 1:
                end = i
            else:
                if end - start >=1 and end!=start:
                    interval.append((start, end))
                start = i
                end = i
        if end - start >=1 and end!=start:
            interval.append((start, end))
    else:
        interval = []
    return interval
