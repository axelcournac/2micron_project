#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 10:54:10 2021
@author: axel
To create a big wig file 
2: change of normalisation and interpolation 
to be synchronised with other algos 
"""
import pyBigWig
import numpy as np
import pandas as pd
import cooler
from scipy import interpolate

# creation of a bw file: 
bw = pyBigWig.open("test.bw", "w")
#bw.addHeader([("chr1", 1000)], maxZooms=0)

bw.addHeader([
("chr1",  230218),
("chr2",  813184),
("chr3",  316620),
("chr4",  1531933),
("chr5",  576874),
("chr6",  270161),
("chr7",  1090940),
("chr8",  562643),
("chr9",  439888),
("chr10",  745751),
("chr11", 666816),
("chr12",  1078177),
("chr13",  924431),
("chr14",  784333),
("chr15",  1091291),
("chr16",  94806)], maxZooms=0)   # keep same order 

sizes = pd.read_table("/home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/sacCer3.chr_sizes.txt",
                      header=None, delimiter="\t") 
zipbObj = zip(sizes[0], sizes[1])
sizes_dist = dict(zipbObj)

list_chr =  ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
             "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
             "chr16"]

# interpolation:
def fill_nan(A):
    '''
    interpolate to fill nan values
    '''
    inds = np.arange(A.shape[0])
    good = np.where(np.isfinite(A))
    f = interpolate.interp1d(inds[good], A[good],bounds_error=False)
    B = np.where(np.isfinite(A),A,f(inds))
    return B

# function for nan sum
def nan_sum(row):
    total = np.nanmean(row)*len(row[~np.isnan(row)])
    return total

# from contact file:
cool_file= "/media/axel/RSG4/diverse_yeast_data/quiescence_2019/fastq/out_Micro-C_WT_log_redone/tmp/valid_idx_pcrfree.pairs.cool"
#cool_file= "/media/axel/RSG4/diverse_yeast_data/AC/fastq_files/seqs/out3_WT_AC_2micron/tmp/valid_idx_pcrfree.pairs.cool"
#cool_file= "/media/axel/RSG4/diverse_yeast_data/FG0089_set2/out_FG0089_Set2/tmp/valid_idx_pcrfree.pairs.2000.cool"


c = cooler.Cooler(cool_file)
cooler.balance_cooler(c, store=True, mad_max=10)   # Normalisation 
name_bank= "Micro-C_WT_log_redone"
BIN=2000
list_chr1 = c.chromnames

chr2 = "plasmid_p2-micron"   

total_number_reads=c.info['sum']
reads_plasmide=0
for chr3 in list_chr1:
            m1 = c.matrix(balance=False).fetch(chr2, chr3)
            if chr2==chr3:
                reads_plasmide += m1.sum()/2+np.trace(m1)/2
            else:
                reads_plasmide += m1.sum()

perc_plasmid=reads_plasmide/total_number_reads
print("Percentage plasmide")
print(perc_plasmid)


for chr1 in list_chr :
    print(chr1)
    matscn = c.matrix(balance=True).fetch(chr1, chr1)
    coverage1 = np.apply_along_axis(nan_sum, 0, matscn)
    mat = c.matrix(balance=True).fetch(chr1, chr2)
    coverage_plasmid1 = np.apply_along_axis(nan_sum, 1, mat)

    # to work with the raw number of reads 
    reads_chr1=0
    coverage=0
    for chr3 in list_chr1:
        m1 = c.matrix(balance=False).fetch(chr1, chr3)
        reads_chr1 += m1.sum()
        coverage_temp=np.apply_along_axis(nan_sum, 1, m1)
        coverage_temp[np.isnan(coverage1)] = np.nan
        coverage += coverage_temp

    m2 = c.matrix(balance=False).fetch(chr2, chr2)
    reads_chr2= m2.sum()
    contact_of_plasmid= np.apply_along_axis(nan_sum, 1, m2)
    m12 = c.matrix(balance=False).fetch(chr1, chr2)
    reads_chr12= m12.sum()
    coverage_plasmid = np.apply_along_axis(nan_sum, 1, m12)
    coverage_plasmid [ np.isnan(coverage_plasmid1)] = np.nan

    coverage_plasmid=coverage_plasmid/coverage # normalisation 
    coverage_plasmid=coverage_plasmid/perc_plasmid # taking into accournt plasmid presence
    
    # interpolation:

    if BIN == 2000:
        proportion_nan = np.count_nonzero(np.isnan(coverage_plasmid))/len(coverage_plasmid)
        if proportion_nan  < 0.9:
            coverage_plasmid_new = fill_nan(coverage_plasmid)
        else :
            coverage_plasmid_new = np.zeros(len(coverage_plasmid))

    
    # writting    
    chroms = np.array([chr1] * len(coverage_plasmid) )
    starts = np.arange(0,len(coverage_plasmid)*BIN,BIN)
    ends = starts+BIN
    values0 = coverage_plasmid

    bw.addEntries(chroms, starts, ends=ends, values=values0)
    
bw.close()





#bw = pyBigWig.open("delete.bw", "w")
#bw.addHeader([("1", 1000)], maxZooms=0)
#chroms = np.array(["1"] * 10)
#starts = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90], dtype=np.int64)
#ends = np.array([5, 15, 25, 35, 45, 55, 65, 75, 85, 95], dtype=np.int64)
#values0 = np.array(np.random.random_sample(10), dtype=np.float64)
#bw.addEntries(chroms, starts, ends=ends, values=values0)
#bw.close()



