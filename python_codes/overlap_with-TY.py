#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 13:33:39 2023
@author: axel
To compute overlab between contacted regions and TY (Transposons of Yeast)
"""

import pandas as pd
import pyranges as pr
import numpy as np
import matplotlib.pyplot as plt


# contact regions by 2mu:
file="/home/axel/Bureau/figure_all_chrms_sup1/redone_norm3/th_08/HSC_plasmids_in_Micro-C_WT_log_SC288_genome.txt.sort2"
df= pd.read_table(file)

#rename of columns:
df.rename(columns = {'chrom1':'Chromosome'}, inplace = True)
df.rename(columns = {'start1':'Start'}, inplace = True)
df.rename(columns = {'end1':'End'}, inplace = True)
        
# TY regions:
file2="/home/axel/Bureau/YEAST/TY/positions_TY.txt"
df2= pd.read_table(file2)

# dictionary of sizes of chromosomes
sizes = pd.read_table("/home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/sacCer3.chr_sizes.txt",
                      header=None, delimiter="\t") 
zipbObj = zip(sizes[0], sizes[1])
sizes_dist = dict(zipbObj)
sizes_dist.popitem()
sizes_dist.popitem()

# conversion into pyranges objects
p1=pr.PyRanges(df)
p2=pr.PyRanges(df2)

len(p1)
len(p2)

# computation of intersections: 
pi=p1.intersect(p2, how="first")
len(pi)

# random group 
sizes = np.array(np.abs(df['Start']-df['End']))
mean_size = int(np.mean(sizes))

list_random=[]
for r in range(1000):
    p_r = pr.random(n=len(p1), length=mean_size, 
                    chromsizes=sizes_dist, int64=False)
    pri=p_r.intersect(p2, how="first")
    list_random.append(len(pri))


len(list_random)

# plot:
plt.hist(list_random,range(np.max(list_random)+2))
plt.xlabel("Number of overlaps")
plt.ylabel("Number of occurences")
plt.title("1000 Random groups")


