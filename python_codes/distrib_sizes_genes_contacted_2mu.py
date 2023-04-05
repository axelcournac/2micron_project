#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 11:21:58 2023
@author: axel
To compare distribution of sizes of genes contacted by 2 mu vs null model
"""

import pandas as pd
import matplotlib.pyplot as plt
import random
import numpy as np
import scipy
 
file="/home/axel/Bureau/figure_all_chrms_sup1/redone_2023/HSC_plasmids_in_Micro-C_WT_log_redone.txt.sort"
df = pd.read_table(file,header=0, delimiter="\t")

plt.hist(df['max_size'],25)
np.mean(df['max_size'])
print(len(df))

genes = pd.read_csv('/home/axel/Bureau/YEAST/GENES_SC288/orf_genomic_all.txt2', 
                    sep=" ", header= None)
size_genes = [ abs(genes[3]-genes[2]) ]
size_genes =size_genes[0]


# random realisation
list_nb_genes=[]
list_max_random=[]
for ri in range(1000):  # N realisations 
    for i in  df['genes']:
        j=i.split(",")
        r=random.choices(size_genes, k=len(j))
        list_nb_genes.append(len(j))
        list_max_random.append(np.max(r))

np.mean(list_nb_genes)

# plot:
plt.hist(df['max_size'],25, label="Contacted area by the 2µ", 
         density=True, fc=(0, 0, 1, 0.5))
plt.hist(list_max_random,25, label="Random null model", 
         density=True, fc=(1, 0, 0, 0.5))

plt.xlabel("Maximum gene size (in bp) in the contacted area (composed of 4 genes in av.)")
plt.ylabel("Density of occurence")

plt.legend()

m1=np.int(np.mean(df['max_size']))
m2=np.int(np.mean(list_max_random))
plt.text(np.mean(df['max_size']), 0.0004 ,"Mean = "+str(m1)+" bp", color=(0, 0, 1, 0.5))
plt.text(np.mean(list_max_random), 0.0003 ,"Mean = "+str(m2)+" bp", color=(1, 0, 0, 0.5))

mann = scipy.stats.mannwhitneyu(x=df['max_size'], y=list_max_random)
plt.text(10000,3,"Mann Whitney p-value=\n"+str(mann.pvalue))
