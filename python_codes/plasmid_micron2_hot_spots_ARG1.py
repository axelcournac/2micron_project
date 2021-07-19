#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 17:08:57 2020
@author: axel
To extract plasmid p2 micron vs chrms from cool files. 
hot spots: to carefully extract hot spots of contact 
between plasmid and host chromosomes 
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
import scn
import ice_mirny3
import distance_law_human

# contact data:
cool_file = sys.argv[1]
name_bank = sys.argv[2]
plasmid_chosen = sys.argv[3]

# cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 /home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/sacCer3.chr_sizes.txt:2000 valid_idx_pcrfree.pairs valid_idx_pcrfree.pairs.cool
c = cooler.Cooler(cool_file)
cooler.balance_cooler(c, store=True, mad_max=5)   # Normalisation 
d=c.info
total_reads = d['sum']
# SC288 assembly:
BIN=2000
list_chr =  ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
             "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
             "chr16","chrM"]

#list_chr = c.chromnames

if len(sys.argv) == 4:
        file_size="/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a91/homologuous_paper/genome_fasta/chrms_sizes.txt2"
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
centro_dist['pYES2'] = 0
centro_dist['pRS416'] = 0
centro_dist['kan'] = 0
centro_dist['laci'] = 0
centro_dist['lacO'] = 0


print(file_size)
print(file_centro)

sizes = pd.read_table(file_size,
                      header=None, delimiter="\t")
zipbObj = zip(sizes[0], sizes[1])
sizes_dist = dict(zipbObj)

#
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
sizes_dist['pST1760'] = 11000
sizes_dist['pYES2'] = 6363
sizes_dist['pRS416'] = 4898
sizes_dist['kan'] = 808
sizes_dist['laci'] = 1083
sizes_dist['lacO'] = 24
#
genes = pd.read_csv('/home/axel/Bureau/YEAST/GENES_SC288/orf_genomic_all.txt2', 
                    sep=" ", header= None)
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
limit_x =  int(1532000/BIN)

if not os.path.exists(name_bank+"_files") :
    os.makedirs(name_bank+"_files")

list_all_contact = []
list_all_chip = []

list_all_contact_lg = []
list_all_chip_lg = []

list_all_chr = []
list_all_pos = []

df_total = pd.DataFrame()

for chr1 in list_chr :
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
    m2 = c.matrix(balance=True).fetch(chr2, chr2)
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
                           height = 0.0006,
                           width = 2.,
                           distance=2)  # for 2000 bp resolution    

    if BIN == 200:
        coverage_plasmid[np.isnan(coverage_plasmid)] = 0.0
        coverage_plasmid_new = coverage_plasmid
        hot_spots = find_peaks(coverage_plasmid_new,
                               height = 0.0011,
                               width = 5.,
                               distance=70)  # for 200 bp resolution

#    plt.plot(coverage_plasmid_new)
#    plt.plot(hot_spots[0],hot_spots[0]*0.0,'o')

    col_chr = [chr1] * len(hot_spots[0])
    col_pos = list(hot_spots[0])
    col_pos = [int(e*BIN) for e in col_pos]
    df = pd.DataFrame({'chrom1': col_chr, 'start1': col_pos})
    df_total = pd.concat([df_total, df])

    # Multiplot:
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(8, 12.8)
    fig.tight_layout(pad=5.0)
    gs = gridspec.GridSpec(4, 1,height_ratios=[7,1,1,1])

    ax1 = plt.subplot(gs[0])
    ax1.imshow(matscn**0.15, interpolation="none", cmap="afmhot_r",
           vmin = 0.0, vmax = 0.8)
    locs, labels = plt.xticks()            # Get locations and labels
    locs2 = [int(i*BIN) for i in locs]
    plt.xticks(locs,locs2)
    fig.tight_layout(pad=2.)
    plt.title("\n\n\n\n"+chr1+" "+str(total_reads)+" "
              +str(reads_chr1)+" "+str(reads_chr2)+" "+str(reads_chr12)+"\n"
              +name_bank)

    ax2 = plt.subplot(gs[1], sharex=ax1)
    ax2.plot(coverage)
    plt.title("Proportion of the signal in INTRA (in %)")

    coverage_plasmid = coverage_plasmid /(float(sizes_dist[chr2]))
    median_plasmid = np.nanmedian(coverage_plasmid)
    print(median_plasmid)
    ax3 = plt.subplot(gs[2], sharex=ax1)
    ax3.plot(coverage_plasmid, color="red" )
    plt.plot(float(centro_dist[chr1])/BIN,0.0,'o', color="orange",label="Centro")
    plt.plot(long_genes_chr[2]/BIN,long_genes_chr[2]*0.0,'>', color="red",label="Long orf")
    plt.plot(long_genes_chr[3]/BIN,long_genes_chr[3]*0.0,'<', color="red",label="Long orf")
#    plot(motif_chr[1]/BIN,motif_chr[1]*0.0,'o', color="red",label="Long orf", alpha=0.1)
#    plot(ori_chr[1]/BIN,ori_chr[1]*0.0,'o', color="yellow",label="Ori")
#    plot(ars_chr[1]/BIN,ars_chr[1]*0.0,'|', color="purple",label="ARS")
#    plot(acetyl_chr[2]/BIN,acetyl_chr[2]*0.0,'|', color="purple",label="acetyl")
    plt.axhline(y=median_plasmid,ls='--',color="Black")
#    plt.ylim(0, limit_y)
#    plt.xlim(0, limit_x)
    plt.plot(hot_spots[0],hot_spots[0]/hot_spots[0]*limit_y,'v', color="teal")
    h=hot_spots[1]
    left_ips = h['left_ips']
    right_ips = h['right_ips']
    plt.plot(left_ips,left_ips*0.0,'|')
    plt.plot(right_ips,right_ips*0.0,'|')
    plt.title(chr2+" Median="+str(round(median_plasmid*10**8,2))+"x 10^8")

    ax4 = plt.subplot(gs[3], sharex=ax1)
    ax4.plot(coverage_plasmid_new, color="green")
#    plt.axhline(y=median_chrM,ls='--',color="Black")
#    plt.ylim(0, limit_y)
#    ax4.plot(b_ip, values_chip, color="royalblue")
    plt.xlabel("Position along the chromosome (bins "+ str(BIN) +"bp)")
#    plt.title(str(chrM)+" Median="+str(round(median_chrM*10**8,2))+"x 10^8")
    plt.title("Plasmid contact with 1D interpolation")

    list_all_contact=np.concatenate((list_all_contact, coverage_plasmid), axis=0)
    list_all_chr=np.concatenate((list_all_chr, [chr1] * len(coverage_plasmid)), axis=0)
    list_all_pos=np.concatenate((list_all_pos, range(len(coverage_plasmid))), axis=0)
    
    # only for long genes:
    l1=long_genes_chr[3]/BIN
    l1=list(l1)
    l1 = [int(x) for x in l1]
    list_all_contact_lg=np.concatenate((list_all_contact_lg,
                                       coverage_plasmid[l1]), axis=0)
    plt.savefig(name_bank+"_files"+"/MAT_SCN_"+chr1+"_"+
            name_bank+"_"+str(BIN/1000)+"kb"+".pdf")

# writting of all hot spots of contact in one file:
df_total.to_csv("HSC_plasmids_in_"+name_bank+".txt", sep='\t',
                index=False)

df_all = pd.DataFrame(
    {'chr': list_all_chr,
     'bin': list_all_pos,
     'contact': list_all_contact
    })

df_all.to_csv("contact_signal_"+chr2+"_"+name_bank+".txt", sep='\t',
                index=False,na_rep='NA')













