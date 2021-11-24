#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 17:08:57 2020
@author: axel
To extract plasmid p2 micron vs chrms from cool files.
3: adding input signal for chip exo 
we remove all part with contact
"""
#import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
sys.path.insert(0, os.path.abspath("/home/axel/Bureau/z_python_scripts_copy"))
import numpy as np
import sys
import re

# SC288 assembly:
BIN=2000
# dictionary :
names = pd.read_table("/media/axel/RSG4/diverse_yeast_data/SraRunTable_ChIP_exo.txt",
                       header=None, delimiter=" ") 
zipbObj = zip(names[0], names[1])
# Create a dictionary from zip object
names_dict = dict(zipbObj)
z=  names[0]
z2 = [i[0:11] for i in z]
z3 = [ names_dict[i] for i in z2]
z4= [re.findall(r"^\d+_(\w+)_\w+_?-?ChIP", i) for i in z3]
zipbObj = zip(names[0], z4)
names2_dict = dict(zipbObj)

# ChIPseq
file_chip = sys.argv[1]
#file_chip='/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR11466946.1_1.fastq.sam.MQ30'
chip_ip = pd.read_csv(file_chip, 
                    sep=" ", header= None)
#chip_input = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR13866613_1_18504_NoTag_BY4741_ChIP-exo.fastq.sam.MQ30', 
#                    sep=" ", header= None)

chip_input = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR13866613.1.1_1.fastq.sam.MQ0', 
                    sep=" ", header= None)

name_prot= sys.argv[1]
name_prot=name_prot[0:11]
name_prot2=names2_dict[name_prot]
name_prot2=name_prot2[0]

list_chr = set(chip_input[0])

def sorted_nicely(l) :
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key)] 
    return sorted(l, key = alphanum_key)

list_chr = sorted_nicely(list_chr)
hist1 = {}
for chr1 in list_chr :       
    chip_ip_chr = chip_ip[chip_ip[0] ==chr1]
    chip_input_chr = chip_input[chip_input[0] ==chr1]
    
    maxi1= max(chip_ip_chr[1]) 
    maxi2= max(chip_input_chr[1]) 
    maxi = max(maxi1, maxi2)
    
    # computation of histograms:
    v_ip, b_ip = np.histogram(chip_ip_chr[1],
                              bins=range(0,maxi+BIN,BIN), density=True)
    v_input, b_input = np.histogram(chip_input_chr[1],
                            bins=range(0,maxi+BIN,BIN), density=True)
    
    b_ip = b_ip/BIN
    b_ip = b_ip[:len(b_ip)-1]
    b_input = b_input/BIN
    b_input = b_input[:len(b_input)-1]
    values_chip = v_ip / v_input
#    values_chip = v_ip  # for RNAseq or Exo ChIP
#    values_chip = np.log(values_chip)  # for RNAseq or other 
    hist1[chr1]=values_chip

#------------------------------------------------------------------------------
# agglomerated 1D plot
#file_set_positions="/home/axel/Bureau/z_python_scripts_copy/Log_2000_files/HSC_plasmids_in_Log_2000.txt"
#file_set_positions="/home/axel/Bureau/YEAST/long_genes.txt2"
#file_set_positions="/home/axel/Bureau/YEAST/GENES_SC288/long_genes_only_host_chrm.txt3"
#file_set_positions="positive_seq.txt2"
file_set_positions="/home/axel/Bureau/figure_all_chrms_sup1/redone_norm3/th_08/HSC_plasmids_in_Micro-C_WT_log_SC288_genome.txt.sort.formated"

pos_set = pd.read_table(file_set_positions)
#name_pos_set="Hot_spots_plasmid"
name_pos_set="73HSC_plasmids_in_Micro-C_WT_log_SC288_genome.txt.sort.formated"

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

print("general score:")
print(score)

# writing of score into a file
fic = open("scores_all_proteins_73HSC_plasmid_ratio_mean_2000bp.txt", "a")
fic.write(name_prot+' '+name_prot2+' '+str(score) +' '+ str(len(chip_ip)) +"\n")
fic.close()

#------------------------------------------------------------------------------
#  Plot of averaged signals :
#s= signal_averaged1 / signal_occ
#plt.plot(s, linewidth = 3.0, label=name_prot,color="red")
#
#plt.plot(signal_mean, linewidth = 3.0, label=name_prot2,color="blue")
#plt.fill_between(range(len(signal_median)),
#                 signal_mean-signal_std/2.0,signal_mean+signal_std/2.0,
#                 alpha=.33,
#                 color = 'blue')
#
#plt.axvline(area, color='k', linestyle='--')
#plt.axhline(1.0, color='grey', linestyle='--')
#
#tick_locs=(0,area,2*area)
#tick_lbls=('-50 kb','Set','+50 kb')
#plt.xticks(tick_locs, tick_lbls,fontsize=10,rotation=45)
#plt.ylabel("ChIP / Input "+name_prot2)
#plt.legend()
#plt.title(str(len(pos_set_uniq))+" "+name_pos_set+"\n"+
#          name_prot2, fontsize=8)
#
#plt.savefig("1D_ENRICHMENT_mean_"+name_prot+"_"+name_prot2+"_"+
#            name_pos_set+".pdf",  dpi=600, format='pdf')
#plt.close("all")



