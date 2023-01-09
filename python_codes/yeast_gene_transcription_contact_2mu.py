#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 17:08:57 2020
@author: axel
To extract plasmid p2 micron contact, transcription level for each gene
"""
import time
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf
import random
import sys
import os
from scipy.stats.stats import spearmanr
sys.path.insert(0, os.path.abspath("/home/axel/Bureau/z_python_scripts_copy"))
import scn 
import ice_mirny3
import scipy
import scipy.ndimage
import scipy.io as sio
import distance_law_human
import hicstuff as hcs
import numpy as np
import json
import sys
import chromosight.utils.plotting as cup
import pandas as pd
import cooler
import re
from Bio import SeqIO
from scipy.stats.stats import pearsonr
import statsmodels.api as sm
lowess = sm.nonparametric.lowess

# contact data:
cool_file= "/media/axel/RSG4/diverse_yeast_data/quiescence_2019/fastq/out_Micro-C_WT_log_redone/tmp/valid_idx_pcrfree.pairs.cool.200"
name_bank= "Log_200"

# cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 /home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/sacCer3.chr_sizes.txt:2000 valid_idx_pcrfree.pairs valid_idx_pcrfree.pairs.cool
c = cooler.Cooler(cool_file)
cooler.balance_cooler(c, store=True, mad_max=10)   # Normalisation 
d=c.info
total_reads = d['sum']
# SC288 assembly:
BIN=200
list_chr =  ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
             "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
             "chr16"]

centro = pd.read_table("/home/axel/Bureau/YEAST/centro1.dat4",
                       header=None, delimiter=" ")
# Create a zip object from two lists
zipbObj = zip(centro[0], centro[1])
# Create a dictionary from zip object
centro_dist = dict(zipbObj)
centro_dist['chrM'] = 0
centro_dist['plasmid_p2-micron'] = 0

sizes = pd.read_table("/home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/sacCer3.chr_sizes.txt",
                      header=None, delimiter="\t") 
zipbObj = zip(sizes[0], sizes[1])
sizes_dist = dict(zipbObj)

#
chr2 = "plasmid_p2-micron"   # plasmid micron
chrM = "chrM"   # mito 

sizes_dist[chrM] = 100000
sizes_dist['plasmid_p2-micron']=6300
sizes_dist['pRS413'] = 4970  # the one of Fabien 
sizes_dist['pSH47'] = 6979
#
genes = pd.read_csv('/home/axel/Bureau/YEAST/GENES_SC288/orf_genomic_all.txt2', 
                    sep=" ", header= None)

genes = genes[genes[0].str.contains("plasmid_p2-micron")==False]
genes = genes[genes[0].str.contains("chrM")==False]

long_genes = genes[ abs(genes[3]-genes[2]) > 7000 ]
len(long_genes)

# ChIPseq
chip_ip = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR7175393.fastq.sam.MQ30.Rpb3_in_Log_All_Replicates_Merged_IP', 
                    sep=" ", header= None)
chip_input = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR7175394.fastq.sam.MQ30.Rpb3_in_Log_All_Replicates_Merged_Input', 
                    sep=" ", header= None)

# Computation of the histo signals 
#-------- 
signal_with_input = 2 # variable telling if there is an input signal 
  
file_set_positions="/home/axel/Bureau/YEAST/GENES_SC288/long_genes_only_host_chrm.txt32"
#file_set_positions="/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/positive_seq.txt2"
#file_set_positions= sys.argv[3]
pos_set = pd.read_table(file_set_positions, sep="\t", header= 0)

name_pos_set="Long_genes"
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
#-------
bin_histo=200
list_chr = set(chip_ip[0])

def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

list_chr = sorted_nicely(list_chr)

#-----
list_all_chip = []  
hist1 = {}
contact_with_plasmid_dict = {}

# Computation of histogram of chip seq signals: 
for chr1 in list_chr : 
    chip_ip_chr = chip_ip[chip_ip[0] ==chr1]
    if signal_with_input != 0 :
        chip_input_chr = chip_input[chip_input[0] ==chr1]
    
    fasta_file="/home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/"+chr1+".fa"
    record = SeqIO.read(open(fasta_file), "fasta") 
    
    maxi1= len(record) 
    if signal_with_input != 0 :
        maxi2= max(chip_input_chr[1]) 
    else :
        maxi2=0    
    maxi = max(maxi1, maxi2)
    
    v_ip, b_ip = np.histogram(chip_ip_chr[1],
                        bins=range(0,maxi+bin_histo,bin_histo), density=True)
    if signal_with_input != 0 :
        v_input, b_input = np.histogram(chip_input_chr[1],
                        bins=range(0,maxi+bin_histo,bin_histo), density=True)
    
    b_ip = b_ip/bin_histo
    b_ip = b_ip[:len(b_ip)-1]
    if signal_with_input != 0 :
        b_input = b_input/bin_histo
        b_input = b_input[:len(b_input)-1]
    
    if signal_with_input != 0 :
        values_chip = v_ip / v_input
    else : 
        values_chip = v_ip  # for RNAseq or Exo ChIP
#        values_chip = np.log(values_chip)  # for RNAseq or other 

    where_are_NaNs = np.isnan(values_chip)
    values_chip[where_are_NaNs] = 0
    list_all_chip=np.concatenate((list_all_chip, values_chip), axis=0)
    hist1[chr1]=values_chip
    
    # contact with plasmid 
    matscn = c.matrix().fetch(chr1, chr2) # ! we focus on the region of contact
    matscn[np.isnan(matscn)] = 0 
    contact_with_plasmid = matscn.sum(axis=1)
    contact_with_plasmid_dict[chr1] = contact_with_plasmid


# Creation of the table gene _ size _ transcription - contact with 2u: 
name_prot="Rpb3"
list_all_contact = []
list_all_chip = []
list_all_sizes = []

genes1=np.array(genes)
for g in range(len(genes1) ) :
    chr1=genes1[g,0]
    pos1=int(genes1[g,2]/BIN)
    pos2=int(genes1[g,3]/BIN)
    if pos1> pos2:
        pos=pos1;pos1=pos2;pos2=pos
        
    signal_gene = hist1[chr1][range(pos1, pos2+1)]
    signal_gene = signal_gene.mean()
        
    contact_plasmid_gene = contact_with_plasmid_dict[chr1][range(pos1, pos2+1)]
#    if len(contact_plasmid_gene>=0) :
    contact_plasmid_gene = np.nanmean(contact_plasmid_gene)
    list_all_contact.append(contact_plasmid_gene)
    list_all_chip.append(signal_gene)
    list_all_sizes.append(abs(genes1[g,3]-genes1[g,2]))
        
    if abs(genes1[g,3]-genes1[g,2]) > 7000:
#        plt.plot(abs(genes1[g,3]-genes1[g,2]),contact_plasmid_gene, 'o', color="orange" )
        plt.text(abs(genes1[g,3]-genes1[g,2]),contact_plasmid_gene, genes1[g,1])
    
#plot   
#    print(genes1[g,1],abs(genes1[g,3]-genes1[g,2]), signal_gene,contact_plasmid_gene) 
        
plt.scatter(list_all_sizes,list_all_contact,
            c=np.log(list_all_chip), cmap=plt.cm.bwr,s =50,
            linewidths=0.5, edgecolors="black", vmin=-1,vmax=1.0)
plt.show()
        
plt.xlabel("Gene size (in bp)")
plt.ylabel("Contact with 2 micron plasmid")

plt.axhline(y=np.mean(list_all_contact),ls='--',color="Black")
plt.colorbar()

PC=pearsonr(list_all_sizes,list_all_contact)

# lowess fir:
z = lowess(list_all_contact,list_all_sizes, 
           frac=1)

plt.plot(z[:,0],z[:,1],lw=3,label="LOWESS FIT", color="black")
plt.legend()

