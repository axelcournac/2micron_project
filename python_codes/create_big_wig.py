#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 10:54:10 2021
@author: axel
To create a big wig file 
"""
import pyBigWig
import numpy as np
import pandas as pd
import cooler

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

## from pairs file
#file_contact="/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a91/quiescence_paper/out_SRR7939017_wt_log/tmp/valid_idx_pcrfree.pairs"
#df=pd.read_table(file_contact,header=None, delimiter="\t", skiprows=21)  # contacts
#
#df2 = df.loc[( (df[3] == "plasmid_p2-micron") & (df[1] !=  "plasmid_p2-micron")) ]
#bin_histo=200
#for chr1 in list_chr :
#    print(chr1)
#    maxi = sizes_dist[chr1]
#    df_chr= df2.loc[df2[1]==chr1]
#    v_ip, b_ip = np.histogram(df_chr[2],
#                              bins=range(0,maxi,bin_histo), 
#                              density=True)
#    b_ip = b_ip[:-1].copy()
#
#    chroms = np.array([chr1] * len(v_ip) )
#    starts = np.array(b_ip, dtype=np.int64)
#    ends = np.array(b_ip+bin_histo, dtype=np.int64)
#    values0 = np.array(v_ip, dtype=np.float64)
#
#    bw.addEntries(chroms, starts, ends=ends, values=values0)
#    
#bw.close()
#
#
#bw2 = pyBigWig.open("/home/axel/test.bw")
#bw2.chroms()
#bw2.header()
#
#

# from contact file:

cool_file= "/media/axel/RSG4/diverse_yeast_data/quiescence_2019/fastq/out_Micro-C_WT_log_redone/tmp/valid_idx_pcrfree.pairs.cool.200"
c = cooler.Cooler(cool_file)
cooler.balance_cooler(c, store=True, mad_max=5)   # Normalisation 
name_bank= "Log_200"
BIN=200

chr2="plasmid_p2-micron"
for chr1 in list_chr :
    print(chr1)
    mat = c.matrix().fetch(chr1, chr2)
    mat[np.isnan(mat)] = 0 
    coverage_plasmid = mat.sum(axis=1)
    coverage_plasmid[coverage_plasmid==0] = np.nan
    
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















