#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 17:08:57 2020
@author: axel
we remove all part with contact
agglo on Mnase-seq signal 
"""
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
sys.path.insert(0, os.path.abspath("/home/axel/Bureau/z_python_scripts_copy"))
import numpy as np
import sys
import re
import itertools

# SC288 assembly:
BIN=2000

# ChIPseq
file_chip = sys.argv[1]
#file_chip="/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/wild-typeSet2-HA_SRR7472982.sam.MQ30"
chip_ip = pd.read_csv(file_chip, 
                    sep=" ", header= None)

name_prot = sys.argv[2]

list_chr =  ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
             "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
             "chr16"]

def sorted_nicely(l) :
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key)] 
    return sorted(l, key = alphanum_key)

#list_chr = sorted_nicely(list_chr)
list_all_values=[]    
hist1 = {}
for chr1 in list_chr :       
    chip_ip_chr = chip_ip[chip_ip[0] ==chr1]
    maxi= max(chip_ip_chr[1]) 
    # computation of histograms:
    v_ip, b_ip = np.histogram(chip_ip_chr[1],
                              bins=range(0,maxi+BIN,BIN), density=True)
    b_ip = b_ip/BIN
    b_ip = b_ip[:len(b_ip)-1]
    values_chip = v_ip 
#    values_chip = v_ip  # for RNAseq or Exo ChIP
    values_chip = np.log(values_chip)  # for RNAseq 
    hist1[chr1]=values_chip
    list_all_values.append(values_chip)

list_all_values= list(itertools.chain.from_iterable(list_all_values))
mean_signal = np.abs(np.mean(np.ma.masked_invalid(list_all_values)))
#------------------------------------------------------------------------------
# agglomerated 1D plot around the spots of contact of 2u plasmid: 
#file_set_positions="/home/axel/Bureau/z_python_scripts_copy/Log_2000_files/HSC_plasmids_in_Log_2000.txt"
#file_set_positions="/home/axel/Bureau/YEAST/long_genes.txt2"
#file_set_positions="/home/axel/Bureau/YEAST/GENES_SC288/long_genes_only_host_chrm.txt3"
file_set_positions="/home/axel/Bureau/figure_all_chrms_sup1/redone_norm3/th_08/HSC_plasmids_in_Micro-C_WT_log_SC288_genome.txt.sort"

pos_set = pd.read_table(file_set_positions)
#name_pos_set="Hot_spots_plasmid"
name_pos_set="HSC_plasmid_in_log_norm3"

input_file_dim = 1   

if input_file_dim == 1 :
   d = {'chrom1': pos_set.chrom1, 'start1': pos_set.start1}
   p1 = pd.DataFrame(data=d)
   p12 = p1
   p12 = p12.drop_duplicates(subset=None, keep='first', inplace=False)

if input_file_dim == 2 :
    d = {'chrom1': pos_set.chrom1, 'start1': pos_set.start1}
    p1 = pd.DataFrame(data=d)
    d = {'chrom1': pos_set.chrom2, 'start1': pos_set.start2}
    p2 = pd.DataFrame(data=d)

    p12  =  pd.concat([p1, p2])
    p12 = p12.drop_duplicates(subset=None, keep='first', inplace=False)

print(len(p12))

len(pos_set)
pos_set_uniq = pos_set.drop_duplicates()  # to remove duplicate entries
len(pos_set_uniq)

area = 20 #  number of bins to compute the signal around the positions
signal_averaged1 = np.zeros(area*2 +1)
signal_averaged2 = np.zeros(area*2 +1)
signal_occ       = np.zeros(area*2 +1)

signal_averaged_list = np.empty((area*2 +1), dtype=np.object)
for i in range(signal_averaged_list.shape[0]):
    signal_averaged_list[i] = []

positions_object=0
for chr1 in list_chr :
    median_chr = np.nanmedian(hist1[chr1])
    print(chr1)
    b = pos_set_uniq.loc[(pos_set_uniq['chrom1'] == chr1)]
    if len(pos_set.columns) > 2:
        b=(b['start1'] + b['end1']) /2
    else: 
        b = b['start1']
    b = uniquearray = list(set(b) )
    for i in range(0, len(b) ) :
        pos_binned =  int(  b[i]  / BIN)
        pi=0
        for p in range(pos_binned-area, pos_binned+area+1):
            if p >=0 and p < len(hist1[chr1]) :
                if np.isfinite(hist1[chr1][p]):
                    signal_averaged1[pi] = signal_averaged1[pi] + hist1[chr1][p]
                    signal_occ[pi] = signal_occ[pi] +1
                    signal_averaged_list[pi].append(hist1[chr1][p])
            pi=pi+1
        positions_object += 1

print(positions_object)

# computation of median and sd signals:
signal_mean = np.zeros(area*2 +1)
signal_median = np.zeros(area*2 +1)
signal_std = np.zeros(area*2 +1)
for i in range(area*2 +1) :
    signal_mean[i] =  np.nanmean(signal_averaged_list[i])
    signal_median[i] =  np.nanmedian(signal_averaged_list[i])
    signal_std[i] =  np.nanstd(signal_averaged_list[i])

# general score renormalised:
score=0
for si in range(area-1,area+2) :
    score = score + signal_mean[si]
score = score / len(range(area-1,area+2))
score=score/mean_signal

print("general score:")
print(score)

# writing of score into a file
fic = open("scores_"+name_prot+"_HSC_plasmid_ratio_mean_2000bp.txt", "a")
fic.write(name_prot+' '+str(score) +' '+ str(len(chip_ip)) +"\n")
fic.close()

#------------------------------------------------------------------------------
#  Plot of averaged signals :
plt.figure(figsize=(3, 2.7))
#s= signal_averaged1 / signal_occ
#plt.plot(s, linewidth = 3.0, label=name_prot,color="red")
#color1="darkblue"
#color1="gold"
color1="blue"
#color1="orange"
plt.plot(signal_mean/mean_signal, linewidth = 3.0, label=name_prot,color=color1)
plt.fill_between(range(len(signal_median)),
                 (signal_mean-signal_std/2.0)/mean_signal,
                 (signal_mean+signal_std/2.0)/mean_signal,
                 alpha=.33,
                 color = color1)

plt.axvline(area, color='k', linestyle='--')
#plt.axhline(1.0, color='grey', linestyle='--')

tick_locs=(0,int(area),int(2*area))
plt.xlim(0,int(2*area))
tick_lbls=('-'+str(area*2)+' kb','Set','+'+str(area*2)+' kb')
plt.xticks(tick_locs, tick_lbls,rotation=0)

plt.ylabel(name_prot)
plt.legend()
plt.title(str(len(pos_set_uniq))+" "+name_pos_set+"\n"+
          name_prot, fontsize=8)
#plt.ylim(0,2.0)
plt.savefig("1D_ENRICHMENT_mean_"+name_prot+"_"+name_prot+"_"+
            name_pos_set+".png",  dpi=600, format='png')

plt.savefig("1D_ENRICHMENT_mean_"+name_prot+"_"+name_prot+"_"+
            name_pos_set+".pdf",  dpi=600, format='pdf')

plt.close("all")


