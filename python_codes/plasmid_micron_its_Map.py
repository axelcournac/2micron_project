#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 17:08:57 2020
@author: axel
To extract plasmid p2 micron vs chrms from cool files.
To generate and annotate the plasmid Map! 
"""
import time
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf
import random
import sys
import os
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
from matplotlib.patches import Rectangle

# %matplotlib qt
number_of_colors = 8

color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]


cool_file="/media/axel/RSG4/diverse_yeast_data/quiescence_2019/fastq/out_Micro-C_WT_log_redone/tmp/valid_idx_pcrfree.pairs.cool.200"
name_bank = "Micro-C_WT_log_redone"

# input data:
cool_file = sys.argv[1]
name_bank = sys.argv[2]
plasmid_chosen = sys.argv[3]

#----------
if not os.path.exists(name_bank+"_files") :
    os.makedirs(name_bank+"_files")

c = cooler.Cooler(cool_file)
cooler.balance_cooler(c, store=True, mad_max=15)   # Normalisation 

 
# SC288 assembly:
BIN=200
centro = pd.read_table("/home/axel/Bureau/YEAST/centro1.dat4",
                       header=None, delimiter=" ") 
# Create a zip object from two lists
zipbObj = zip(centro[0], centro[1])
# Create a dictionary from zip object
centro_dist = dict(zipbObj)
centro_dist['chrM'] = 0
centro_dist['plasmid_p2-micron'] = 0

sizes = pd.read_table("/home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron_pYADE4/sacCer3.chr_sizes.txt",
                      header=None, delimiter="\t") 
zipbObj = zip(sizes[0], sizes[1])
sizes_dist = dict(zipbObj)

#
#chr2 = "plasmid_p2-micron"   # plasmid micron
if plasmid_chosen == "None":
    chr2 = "plasmid_p2-micron"   
else : 
    chr2 = plasmid_chosen

list_chr =  [plasmid_chosen]

chrM = "chrM"   # mito 
chrM = ('chr2', 250000, 350000) 
sizes_dist[chrM] = 100000
sizes_dist['pST1760'] = 11000
#
genes = pd.read_csv('/home/axel/Bureau/YEAST/GENES_SC288/orf_genomic_all.txt2', 
                    sep=" ", header= None)
long_genes = genes[ abs(genes[3]-genes[2]) > 7000 ]
len(long_genes)

plasmid_genes = genes[ genes[0] == "plasmid_p2-micron" ]
plasmid_genes = pd.read_csv('/home/axel/Bureau/YEAST/GENES_SC288/genes_plasmid_J01347.1.txt.recorrected', 
                    sep=" ", header= None)

# Genomic data
ori = pd.read_csv('/home/axel/Bureau/YEAST/origins/origins_Alvino.txt2', 
                    sep=" ", header= None)
ars = pd.read_csv('/home/axel/Bureau/YEAST/ARS/ars.txt2', 
                    sep=" ", header= None)
acetyl = pd.read_csv('/home/axel/Bureau/YEAST/H4_acetylation/GSM3445778_YSB787-1-H4ac_peaks.narrowPeak3', 
                    sep=" ", header= None)
motif = pd.read_csv('/home/axel/Bureau/2micron_plasmid_PROJECT/TGCATTTTT.sam2', 
                    sep=" ", header= None)

csfont = {'fontname':'Times New Roman'}
hfont = {'fontname':'Times New Roman'}

limit_y =  8*10**-7
for chr1 in list_chr :
    long_genes_chr = long_genes[long_genes[0] ==chr1]
    ori_chr = ori[ori[0] ==chr1]
    ars_chr = ars[ars[0] ==chr1]
    acetyl_chr = acetyl[acetyl[0] ==chr1]
    
    motif_chr = motif[motif[0] =="plasmid_p2-micron"]
    v, b = np.histogram(motif_chr[1],bins=range(0,max(motif_chr[1]),BIN) )
    b = b/BIN
    b = b[:len(b)-1]
    
    matscn = c.matrix().fetch(chr1, chr1) 
    matscn[np.isnan(matscn)] = 0 
    coverage = matscn.sum(axis=0)
    coverage[coverage==0] = np.nan

    mat = c.matrix().fetch(chr1, chr2)
    mat[np.isnan(mat)] = 0 
    coverage_plasmid = mat.sum(axis=1)
    coverage_plasmid[coverage_plasmid==0] = np.nan
    
    mat2 = c.matrix().fetch(chr1, chrM)
    mat2[np.isnan(mat2)] = 0   
    coverage_chrM = mat2.sum(axis=1)
    coverage_chrM[ coverage_chrM==0 ] = np.nan
    
    # Multiplot: 
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(8, 11.8)
    fig.tight_layout(pad=2.5)
    gs = gridspec.GridSpec(4, 1,height_ratios=[9,1,1,1])

    ax1 = plt.subplot(gs[0])
#    ax1.imshow(matscn**0.5,interpolation="none", cmap="afmhot_r")
    ax1.imshow(matscn**0.4,interpolation="none", cmap="hot_r")
    plt.title(chr1)
    
    ax2 = plt.subplot(gs[1], sharex=ax1)
    ax2.plot(coverage)
    plt.title("Proportion of the signal in INTRA (in %)")
    
    coverage_plasmid = coverage_plasmid /(sizes_dist[chr2])
    median_plasmid = np.nanmedian(coverage_plasmid) 
    print(median_plasmid)
    ax3 = plt.subplot(gs[2], sharex=ax1)
    plt.ylim(-0.25, 0.25)
    plt.title("Genes")
#    plt.title(chr2+" Median="+str(round(median_plasmid*10**8,2))+"x 10^8")
    
    name_gene=np.array(plasmid_genes[1])
    start=np.array(plasmid_genes[2])
    end=np.array(plasmid_genes[3])
    
    for i in range(len(start)-3):
        color_chosen=color[i]
        if name_gene[i] == "Flp":
            color_chosen = "royalblue"
        rectangle = plt.Rectangle(   (min(start[i],end[i])/BIN,0.0),
                                  abs(end[i]-start[i])/BIN, 0.2 , 
                                  fc=color_chosen)
        plt.gca().add_patch(rectangle)
        plt.text(min(start[i],end[i])/BIN, y=-0.1, s=name_gene[i])
        
    for i in range(len(start)-3,len(start)-1):
        rectangle = plt.Rectangle(   (min(start[i],end[i])/BIN,-0.17),
                                  abs(end[i]-start[i])/BIN, 0.1, 
                                  fc='yellow')
        plt.gca().add_patch(rectangle)
        plt.text(min(start[i],end[i])/BIN, y=-0.22, s=name_gene[i])    
        
#    plt.tick_params(
#    axis='x',          # changes apply to the x-axis
#    which='both',      # both major and minor ticks are affected
#    bottom=False,      # ticks along the bottom edge are off
#    top=False,         # ticks along the top edge are off
#    labelbottom=False) # labels along the bottom edge are off
        
    coverage_chrM = coverage_chrM/(sizes_dist[chrM]) 
    median_chrM = np.nanmedian(coverage_chrM) 
    ax4 = plt.subplot(gs[3], sharex=ax1)
#    ax4.plot(coverage_chrM, color="green")
    ax4.plot(b, v ,'o-',color="royalblue")
#    plt.axhline(y=median_chrM,ls='--',color="Black")
#    ax4.plot(b, v, color="royalblue")
#    plt.xlabel("Position along the chromosome (bins 200 bp)")
#    plt.title(str(chrM)+" Median="+str(round(median_chrM*10**8,2))+"x 10^8")
#    plt.title("Genes")
    plt.xlabel("Position along the chromosome (bins 200 bp)")
    plt.title("TGCATTTTT density")
    locs, labels = plt.xticks()            # Get locations and labels
    locs2 = [int(i*BIN) for i in locs]
    plt.xticks(locs,locs2)
    fig.tight_layout(pad=2.)
    plt.savefig(name_bank+"_files"+'/MAT_SCN_'+chr1+"_"+
            name_bank+"_"+str(BIN/1000)+"kb"+'.pdf')
    
    
    
