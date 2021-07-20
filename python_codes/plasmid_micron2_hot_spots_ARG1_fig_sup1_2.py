#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 17:08:57 2020
@author: axel
To extract plasmid p2 micron vs chrms from cool files. 
hot spots: to carefully extract hot spots of contact 
between plasmid and host chromosomes 
to carry out in spyder
"""
import time
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf
import random as rando
import sys
import os

import scipy
import scipy.ndimage
import scipy.io as sio
from scipy.signal import find_peaks
from scipy import interpolate

import numpy as np
import json
import cooler

import chromosight.utils.plotting as cup 
import hicstuff as hcs
sys.path.insert(0, os.path.abspath("/home/axel/Bureau/z_python_scripts_copy"))


# contact data:
cool_file = sys.argv[1]
name_bank = sys.argv[2]
plasmid_chosen = sys.argv[3]

cool_file="/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a91/quiescence_paper/fastq_sam/out_Micro-C_WT_log_SC288_genome/tmp/valid_idx_pcrfree.pairs.cool.2000"
name_bank = "Micro-C_WT_log_SC288_genome"
plasmid_chosen = 'plasmid_p2-micron'

# cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 /home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/sacCer3.chr_sizes.txt:2000 valid_idx_pcrfree.pairs valid_idx_pcrfree.pairs.cool
c = cooler.Cooler(cool_file)
cooler.balance_cooler(c, store=True, mad_max=10)   # Normalisation 
d=c.info
total_reads = d['sum']
# SC288 assembly:
BIN=2000
list_chr =  ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
             "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
             "chr16"]

list_chr = c.chromnames

if len(sys.argv) == 4:
        file_size="/home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron_pYADE4/sacCer3.chr_sizes.txt2"
        file_centro="/home/axel/Bureau/YEAST/centro1.dat5"
else :
        file_size = sys.argv[4]
        file_centro = sys.argv[5]

centro = pd.read_table(file_centro,
                       header=None, delimiter="\t")
# Create a zip object from two lists
zipbObj = zip(centro[0], centro[1])
# Create a dictionary from zip object
centro_dist = dict(zipbObj)
centro_dist['chrM'] = 0
centro_dist['plasmid_p2-micron'] = 0
centro_dist['pJL638'] = 0

print(file_size)
print(file_centro)

sizes = pd.read_table(file_size,
                      header=None, delimiter="\t")
zipbObj = zip(sizes[0], sizes[1])
sizes_dist = dict(zipbObj)

#
list_chr =  ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
             "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
             "chr16"]

chr2 = "plasmid_p2-micron"   # plasmid micron
chr2="pJL638"
if plasmid_chosen == "None":
    chr2 = "plasmid_p2-micron"
else :
    chr2 = plasmid_chosen

print("Signal of contact computed with the following molecule:")
print(chr2)
chrM = "chrM"   # mito

sizes_dist['pRS413'] = 4970  # the one of Fabien 
sizes_dist['pSH47'] = 6979
sizes_dist['pJL638'] = 7500
#
genes = pd.read_csv('/home/axel/Bureau/YEAST/GENES_SC288/orf_genomic_all.txt2', 
                    sep=" ", header= None)

#genes = pd.read_csv('/home/axel/Bureau/YEAST/GENES_SC288/orf_genomic_all.fasta.infos.reformated.renamed2', 
#                    sep=" ", header= None)

long_genes = genes[ abs(genes[3]-genes[2]) > 7000 ]
len(long_genes)

# Genomic data
ori = pd.read_csv('/home/axel/Bureau/YEAST/origins/origins_Alvino.txt2', 
                    sep=" ", header= None)
ars = pd.read_csv('/home/axel/Bureau/YEAST/ARS/ars.txt2', 
                    sep=" ", header= None)
acetyl = pd.read_csv('/home/axel/Bureau/YEAST/H4_acetylation/GSM3445778_YSB787-1-H4ac_peaks.narrowPeak3', 
                    sep=" ", header= None)
motif = pd.read_csv('/home/axel/Bureau/2micron_plasmid_PROJECT/TGCATTTTT.sam2', 
                    sep=" ", header= None)
#
limit_y =  3.0*10**-7
limit_y =  0.0023
limit_x =  int(1532000/BIN)

if not os.path.exists(name_bank+"_files") :
    os.makedirs(name_bank+"_files")

list_all_contact = []
list_all_chip = []

list_all_contact_lg = []
list_all_chip_lg = []

df_total = pd.DataFrame()

fig2 = matplotlib.pyplot.gcf()
fig2.set_size_inches(8, 12.8)
fig2.tight_layout(pad=5.0)
gs2 = gridspec.GridSpec(17, 1)

fig2.patch.set_visible(False)

jj=0
for chr1 in list_chr :
    genes_chr = genes[genes[0] ==chr1]
    long_genes_chr = long_genes[long_genes[0] ==chr1]
    ori_chr = ori[ori[0] ==chr1]
    ars_chr = ars[ars[0] ==chr1]
    acetyl_chr = acetyl[acetyl[0] ==chr1]

    matscn = c.matrix().fetch(chr1, chr1)
    matscn[np.isnan(matscn)] = 0
    coverage = matscn.sum(axis=0)
    coverage[coverage==0] = np.nan

    mat = c.matrix().fetch(chr1, chr2)
    mat[np.isnan(mat)] = 0
    coverage_plasmid = mat.sum(axis=1)
    coverage_plasmid[coverage_plasmid==0] = np.nan

    m1 = c.matrix(balance=False).fetch(chr1, chr1)
    reads_chr1= m1.sum()
    m2 = c.matrix(balance=False).fetch(chr2, chr2)
    reads_chr2= m2.sum()
    m12 = c.matrix(balance=False).fetch(chr1, chr2)
    reads_chr12= m12.sum()

    # computation of correlation:
    maxi= matscn.shape[0]*BIN

    # !!  Identification of peaks for hot spots of contact:
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

    if BIN == 2000:
        proportion_nan = np.count_nonzero(np.isnan(coverage_plasmid))/ len(coverage_plasmid)
        if proportion_nan  < 0.9:
            coverage_plasmid_new = fill_nan(coverage_plasmid)
        else :
            coverage_plasmid_new = np.zeros(len(coverage_plasmid))
        hot_spots = find_peaks(coverage_plasmid_new, 
                               height = 0.0010,
                               width = 2.,
                               distance=2)  # for 2000 bp resolution

    if BIN == 200:
        coverage_plasmid[np.isnan(coverage_plasmid)] = 0.0
        coverage_plasmid_new = coverage_plasmid
        hot_spots = find_peaks(coverage_plasmid_new,
                               height = 0.0011,
                               width = 5.,
                               distance=70)  # for 200 bp resolution

    col_chr = [chr1] * len(hot_spots[0])
    col_pos = list(hot_spots[0])
    col_pos = [int(e*BIN) for e in col_pos]
    
    col_left= hot_spots[1]["left_ips"]
    col_right= hot_spots[1]["right_ips"]
    col_left = [int(e*BIN) for e in col_left]
    col_right = [int(e*BIN) for e in col_right]
    col_sizes = [e2-e1 for e1,e2 in zip(col_left,col_right)]
    
    heights = hot_spots[1]['peak_heights']
 
    genes_present = [genes_chr.loc[(( ((genes_chr[2]+genes_chr[3])/2.0>e1) & 
                                     ((genes_chr[2]+genes_chr[3])/2.0<e2)) |
                     ( (genes_chr[3]>e1) & (genes_chr[3]<e2))), [1] ] 
    for e1,e2 in zip(col_left,col_right)]
    
    pos1_genes_present = [genes_chr.loc[(( ((genes_chr[2]+genes_chr[3])/2.0>e1) & 
                                     ((genes_chr[2]+genes_chr[3])/2.0<e2)) |
                     ( (genes_chr[3]>e1) & (genes_chr[3]<e2))), [3] ] 
    for e1,e2 in zip(col_left,col_right)]
    
#    pos2_genes_present = [genes_chr.loc[(( (genes_chr[2]>e1) & (genes_chr[2]<e2)) |
#                     ( (genes_chr[3]>e1) & (genes_chr[3]<e2))), [2] ]  
#    for e1,e2 in zip(col_left,col_right)]
    
    pos2_genes_present = [genes_chr.loc[(( ((genes_chr[2]+genes_chr[3])/2.0>e1) & 
                                     ((genes_chr[2]+genes_chr[3])/2.0<e2)) |
                     ( (genes_chr[3]>e1) & (genes_chr[3]<e2))), [2] ] 
    for e1,e2 in zip(col_left,col_right)]
    
      
    sizes_genes_present = [ abs(e1[2]-e2[3]) for e1,e2 in 
                          zip(pos2_genes_present,pos1_genes_present)]
    
    genes_present = [ list(e[1]) for e in genes_present]
    sizes_genes_present = [ list(e) for e in sizes_genes_present]

    
#    df = pd.DataFrame({'chrom1': col_chr, 'start1': col_pos})
    df = pd.DataFrame({'chrom1': col_chr, 'start1': col_left, 
                       'end1': col_right,'size': col_sizes,
                       'heights': heights,
                       'genes': genes_present,
                       'sizes':sizes_genes_present})
    df_total = pd.concat([df_total, df])
    
    coverage_plasmid = coverage_plasmid /(float(sizes_dist[chr2]))
    median_plasmid = np.nanmedian(coverage_plasmid_new)
     
    list_all_contact=np.concatenate((list_all_contact, coverage_plasmid_new), axis=0)
#
    aax1 = plt.subplot(gs2[jj])

    aax1.plot(coverage_plasmid_new, color="darkblue")
    plt.plot(float(centro_dist[chr1])/BIN,0.0,'o', color="red", markersize=4)
    
    # plot with rectangles: 
    pos1=np.array(long_genes_chr[2])
    pos2=np.array(long_genes_chr[3])
    for i in range(len(pos1)):
        rectangle = plt.Rectangle( (min(pos1[i],pos2[i])/BIN,-0.001),
                                  abs(pos2[i]-pos1[i])/BIN, 
                                  0.0011, 
                                  fc="grey")
        plt.gca().add_patch(rectangle)
    
    
    if len(long_genes_chr[0])>0:
        x=np.array(long_genes_chr[3]/BIN)
        y=np.array(long_genes_chr[3]*0.0)
        title=np.array(long_genes_chr[1])
        for z in range(len(long_genes_chr[0])) :
            plt.text(x[z],-0.0012,title[z])
    plt.axhline(y=median_plasmid,xmin=0, xmax=(sizes_dist[chr1]/BIN)/limit_x,
                ls='--',color="Black")
    
    aax1.spines['right'].set_visible(False)
    aax1.spines['top'].set_visible(False)
    aax1.spines['bottom'].set_visible(False)
    
    plt.axvline(len(coverage_plasmid_new),color="Black", lw=1)

    plt.axhline(y=-limit_y*0.01,xmin=0,xmax=(sizes_dist[chr1]/BIN)/limit_x, color="black", lw=1)
    plt.axhline(limit_y-limit_y*0.01, xmin=0,xmax=(sizes_dist[chr1]/BIN)/limit_x, color="black", lw=1)
    
    plt.ylim(-limit_y*0.1, limit_y)
    plt.xlim(0, limit_x)
#    aax1.set_yticks([], [])
    aax1.set_xticks([], [])
    plt.ylabel(chr1)
    plt.plot(hot_spots[0],hot_spots[0]/hot_spots[0]*limit_y,'v', color="orange",
             markersize=10.0)
    h=hot_spots[1]
    left_ips = h['left_ips']
    right_ips = h['right_ips']
#    plt.plot(left_ips,left_ips*0.0,'|')
#    plt.plot(right_ips,right_ips*0.0,'|')
#    plt.title(chr2+" Median="+str(round(median_plasmid*10**8,2))+"x 10^8")
    jj=jj+1

plt.plot(float(centro_dist[chr1])/BIN,0.0,'o', color="red", markersize=4,
         label="Centromere")
plt.plot(hot_spots[0],hot_spots[0]/hot_spots[0]*limit_y,'v', color="orange",
             markersize=10.0, label="Automatically detected")
plt.legend()     

plt.savefig(name_bank+"_files"+"/all_chrs_"+"_"+
        name_bank+"_"+str(BIN/1000)+"kb"+".pdf")



# writting of all hot spots of contact in one file:
df_total_sorted = df_total.sort_values(by=['heights'], ascending=False)

df_total_sorted.to_csv("HSC_plasmids_in_"+name_bank+".txt.sort", sep='\t',
                index=False)

list_all_contact=pd.DataFrame(list_all_contact)

list_all_contact.to_csv("contact_signal_"+chr2+"_"+name_bank+".txt", sep='\t',
                index=False,na_rep='NA')

