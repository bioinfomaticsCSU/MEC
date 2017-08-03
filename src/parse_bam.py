#!/usr/bin/env python

import pysam
import collections
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

def GCProportion(line):
    c = line.count("g") + line.count("c")+line.count("G") + line.count("C")
    t = len(line)*1.0
    if t==0: return 0
    return c/t

def get_all_interval(samfile,minq,alpha,beta,gamma,seq,GCavgRate):

    
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
    print s,d
    print alpha,beta,gamma
    interval = {}
    
    for k,v in zip(references,lengths):
        reads = samfile.fetch(k, 0, v)
        interval_valid = get_valid_interval(k,v,lengths,reads,s,d,minq,alpha,beta,gamma,seq,GCavgRate)
        interval[k] = interval_valid

    return interval

def get_valid_interval(cid,length,lengths,reads,mu,sigma,minq,alpha,beta,gamma,seq,GCavgRate):

    left_accordant = []
    left_discordant = []
    right_discordant = []
    frag_coverage = [0]*length
    
    for read in reads:
        if read.cigarstring!=None:
            if read.rnext == read.tid:
                if mu<1000:
                    if read.is_proper_pair:
                        if read.next_reference_start < read.reference_start:
                            left_accordant.append((read.reference_name,read.next_reference_start,read.reference_name,read.reference_end,read.isize))
                            for i in xrange(read.next_reference_start,read.reference_end):
                                frag_coverage[i]+=1
                else:
                    if mu-3*sigma<= abs(read.isize) <= mu+3*sigma:
                        if read.next_reference_start < read.reference_start:
                            left_accordant.append((read.reference_name,read.next_reference_start,read.reference_name,read.reference_end,read.isize))
                            for i in xrange(read.next_reference_start,read.reference_end):
                                frag_coverage[i]+=1
            else:
                if read.is_reverse:
                    right_discordant.append((read.next_reference_name,read.pnext,read.reference_name,read.pos,read.isize))
                else:
                    left_discordant.append((read.reference_name,read.pos,read.next_reference_name,read.pnext,read.isize))
                                 
    frag_threshold = sum(frag_coverage)*1.0/length
    frag_pos = [i for i in xrange(length) if frag_coverage[i]>alpha*frag_threshold]
    frag_interval = pos_to_interval(frag_pos,100)

    interval = []
    if len(frag_interval)>=2:
        interval_len = len(frag_interval)
        p0 = 0
        for i in xrange(interval_len):
            if i==interval_len-1:
                p1 = length-1
                interval.append((p0,p1))
            else:
                error_pos = [frag_coverage[j] for j in xrange(frag_interval[i][1]+1,frag_interval[i+1][0])]
                break_pos01 = error_pos.index(min(error_pos))
                p1 = frag_interval[i][1]+break_pos01
                interval.append((p0,p1))
                p0 = p1+1
                p1 = frag_interval[i+1][1]
    else:
        interval = [(0,length-1)]
        
    interval_valid=[]
    discordant=[]
    after_discordant=[]
    before_discordant=[]
    accordant=[]
    discordante_rate=[0]*(len(interval)-1)
    accordante_rate=[0]*(len(interval)-1)
    if len(interval)>=2:
        p0 = interval[0][0]
        p1 = interval[0][1]
        interval_len = len(interval)
        for i in xrange(interval_len-1):
            before_ctg =  []
            after_ctg = []
            
            start = max(interval[i][1]-mu,interval[i][0])
            end = min(interval[i+1][0]+mu,interval[i+1][1])
            
            after_discordant.append([item for item in left_discordant if (item[1]>=start and item[1]<interval[i+1][0]) ])
            before_discordant.append([item for item in right_discordant if (item[3]>=interval[i][1]+1 and item[3]<end) ])
            accordant.append([item for item in left_accordant if (item[1]<interval[i][1] and item[3]>interval[i][1]+1) ])

            for item in left_discordant:
                if (item[1]>=start and item[1]<interval[i+1][0]):
                    after_ctg.append(item[2])

            for item in right_discordant:
                if (item[3]>=interval[i][1]+1 and item[3]<end):
                    before_ctg.append(item[0])

            couter01 = collections.Counter(after_ctg)
            couter02 = collections.Counter(before_ctg)
            
            n1=n2=0
            if len(couter01.most_common(1))!=0:
                n1 = couter01.most_common(1)[0][1]
            if len(couter02.most_common(1))!=0:
                n2 = couter02.most_common(1)[0][1]

            isize_discordant = 0
            if n1>=n2:
                isize_discordant += n1
                discordant.append(after_discordant[i])
            else:
                isize_discordant += n2
                discordant.append(before_discordant[i])
            
            isize_accordant = len(accordant[i])                    
            accordante_rate[i] = isize_accordant*1.0/max((isize_accordant+isize_discordant),0.01)
            discordante_rate[i] = isize_discordant*1.0/max((isize_accordant+isize_discordant),0.01)
                
            if accordante_rate[i] >= beta or discordante_rate[i]==0.0:
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
        interval_valid = interval

    intervals=[]
    if len(interval_valid)>=2:
        interval_len = len(interval_valid)
        p0 = interval_valid[0][0]
        p1 = interval_valid[0][1]
        for i in xrange(interval_len-1):
            start = max(interval_valid[i][1]-mu/2,0)
            end = min(interval_valid[i+1][0]+mu/2,length-1)
            
            p_gc = GCProportion(seq[cid][start:end])
            l_gc = GCProportion(seq[cid][start:interval_valid[i][1]])
            r_gc = GCProportion(seq[cid][interval_valid[i+1][0]+1:end])
            
            min_gc = GCavgRate-gamma
            max_gc = GCavgRate+gamma
            if (min(p_gc,l_gc,r_gc)<=min_gc or max(p_gc,l_gc,r_gc)>=max_gc):
                pos = i+1
                p1 = interval_valid[pos][1]
                if pos == interval_len-1:
                    intervals.append((p0,p1))
            else:
                intervals.append((p0,p1))
                pos= i+1
                p0 = interval_valid[pos][0]
                p1 = interval_valid[pos][1]
                if pos == interval_len-1:
                    intervals.append((p0,p1))
    else:
        intervals=interval_valid
            
    return intervals
                
def pos_to_interval(pos_list,minlen):
    
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
                if end - start >minlen:
                    interval.append((start, end))
                start = i
                end = i
        if end - start >minlen:
            interval.append((start, end))
    else:
        interval = []
    return interval
