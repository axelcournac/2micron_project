#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 17:08:57 2020
@author: axel
to extract fasta sequences for loops   
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
from Bio import SeqIO
from Bio.SeqUtils import GC
import domainogram


# loop file from Chromosight 
loop_file = "/home/axel/Bureau/analyse_sequence_HSC/HSC_plasmids_in_SRR7939017_wt_log.txt2" 
loop_file = "/home/axel/Bureau/figure_all_chrms_sup1/redone_norm3/th_08/HSC_plasmids_in_Micro-C_WT_log_SC288_genome.txt.sort.formated2"
group2=pd.read_table(loop_file, delimiter="\t") 
print(len(group2))

fh1 = open("sequences_73HSC.txt","w")

list_loop = []
for index, row in group2.iterrows():
    chr1= row['chrom1']
    start1 =  int(row['start1'])
    end1=  int(row['end1'] )  
    ele1 = chr1+" "+str(start1)+" "+str(end1)

    fasta_file="/home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/"+chr1+".fa"
    record = SeqIO.read(open(fasta_file), "fasta")
    if start1 < end1:
        sequence = record.seq[ start1 : end1]
    else: 
        sequence = record.seq[ end1 : start1]
    fh1.write(">HSC_"+str(chr1)+"_"+ str(start1) + "_" + str(end1) + "\n" )
    fh1.write(str(sequence)+"\n")
         
fh1.close()































 
# Loops or hot spots of contact of plasmid: 
pos_set = pd.read_table('/home/axel/Bureau/YEAST/strong_loops_cohesin_quant_cohe_classic_10-50kb_49_mitotic.txt')


# SC288 assembly:
BIN=200
list_chr =  ["plasmid_p2-micron"]

list_chr =  ["chr12"]

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
area=60000
chr1 = ('chr11', 535647-area, 547925+area) 
list_chr = [chr1]

#chr2 = ('chr12',   306855-area, 315731+area)
#chr2 = ('chr12',   349006-area, 363738+area)
#chr2 = ('chr12',  1043998-area,  1051884+area)


chr2 = "plasmid_p2-micron"   # plasmid micron
#chrM = "chrM"   # mito 

chrM = 'chrM' 

sizes_dist[chrM] = area 
sizes_dist['plasmid_p2-micron']=6300
sizes_dist['pRS413'] = 4970  # the one of Fabien 
sizes_dist['pSH47'] = 6979

#
genes = pd.read_csv('/home/axel/Bureau/YEAST/GENES_SC288/orf_genomic_all.txt2', 
                    sep=" ", header= None)
long_genes = genes[ abs(genes[3]-genes[2]) > 7000 ]
len(long_genes)

plasmid_genes = genes[ genes[0] == "plasmid_p2-micron" ]
plasmid_genes = pd.read_csv('/home/axel/Bureau/YEAST/GENES_SC288/genes_plasmid_J01347.1.txt', 
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

#---------
list_all_contact = np.zeros(int(sizes_dist[chr1]/BIN) +1)
BIN=2000
limit_y =  8*10**-7
for i in range(len(p12)) :
    print(i)
    chr1 = p12.chrom1[i]
    fasta_file="/home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/"+chr1+".fa"
    record = SeqIO.read(open(fasta_file), "fasta")  
    
    pos1 = p12.start1[i]-area
    pos2 = p12.start1[i]+area
    if pos1<0:
        pos1=0
    if pos2>=len(record):
        pos2=len(record)        
    
    long_genes_chr = long_genes[long_genes[0] ==chr1]
    ori_chr = ori[ori[0] ==chr1]
    ars_chr = ars[ars[0] ==chr1]
    acetyl_chr = acetyl[acetyl[0] ==chr1]
    
    motif_chr = motif[motif[0] ==chr1]
    
    v, b = np.histogram(motif_chr[1],bins=range(0,max(motif_chr[1]),BIN) )
    b = b/BIN
    b = b[:len(b)-1]
    
    chr1_focus = (chr1,pos1,pos2) # region of focus  
    
    matscn = c.matrix().fetch(chr1_focus, chr1_focus) # we focus on the region of contact
    matscn[np.isnan(matscn)] = 0 
    coverage = matscn.sum(axis=0)
    coverage[coverage==0] = np.nan

    mat = c.matrix().fetch(chr1_focus, chr2)
    mat[np.isnan(mat)] = 0 
    coverage_plasmid = mat.sum(axis=1)
    coverage_plasmid[coverage_plasmid==0] = np.nan
    
    mat2 = c.matrix().fetch(chr1_focus, chr2)
    mat2[np.isnan(mat2)] = 0   
    coverage_chrM = mat2.sum(axis=1)
    coverage_chrM[coverage_chrM==0] = np.nan
    
    #---------------------
    # Multiplot: 
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(12, 11.8)
    fig.tight_layout(pad=2.5)
    gs = gridspec.GridSpec(2, 2,height_ratios=[9,1],
                           width_ratios=[9,2])

    ax1 = plt.subplot(gs[0])
    ax1.imshow(matscn**0.15,interpolation="none", cmap="afmhot_r", 
             aspect="auto")
    plt.ylabel(chr1+' '+str(int(p12.start1[i])))
    plt.xlabel(chr2)
    plt.title(chr1)
    
    ax11 = plt.subplot(gs[1],sharey=ax1)  # vertical
    v1 = range(mat.shape[0])
    v2 = mat.sum(axis=1)
    plt.title("Contacts with \n"+chr2)
    ax11.plot(v2,v1, linewidth=3.0)
    plt.xticks(rotation=70)
    
    ax2 = plt.subplot(gs[2], sharex=ax1)
    plt.title("Contact with plasmid "+chr2) 
    coverage_plasmid = coverage_plasmid /(sizes_dist[chr2])
    median_plasmid = np.nanmedian(coverage_plasmid) 
    ax2.plot(coverage_plasmid, linewidth=3.0)
    print(median_plasmid)
    
    plt.savefig(name_bank+"_files"+'_'+chr1+'_'+str(pos1)+
                "_"+name_bank+"_"+str(BIN/1000)+"kb_intra"+'.pdf')
    
    
#    ax3 = plt.subplot(gs[4], sharex=ax1)
#    plt.ylim(-0.25, 0.25)
#    plt.title("Plasmid Genes")
##    plt.title(chr2+" Median="+str(round(median_plasmid*10**8,2))+"x 10^8")
#    
#    name_gene=np.array(plasmid_genes[1])
#    start=np.array(plasmid_genes[2])
#    end=np.array(plasmid_genes[3])
#    
#    for i in range(len(start)-3):
#        color_chosen=color[i]
#        if name_gene[i] == "Flp":
#            color_chosen = "royalblue"
#        rectangle = plt.Rectangle(   (min(start[i],end[i])/BIN,0.0),
#                                  abs(end[i]-start[i])/BIN, 0.2 , 
#                                  fc=color_chosen)
#        plt.gca().add_patch(rectangle)
#        plt.text(min(start[i],end[i])/BIN, y=-0.1, s=name_gene[i])
#        
#    for i in range(len(start)-3,len(start)-1):
#        rectangle = plt.Rectangle(   (min(start[i],end[i])/BIN,-0.17),
#                                  abs(end[i]-start[i])/BIN, 0.1, 
#                                  fc='yellow')
#        plt.gca().add_patch(rectangle)
#        plt.text(min(start[i],end[i])/BIN, y=-0.22, s=name_gene[i])    
#        
#  
#    coverage_chrM = coverage_chrM/(sizes_dist[chrM]) 
#    median_chrM = np.nanmedian(coverage_chrM) 
#    
#    ax4 = plt.subplot(gs[6], sharex=ax1)
##    ax4.plot(coverage_chrM, color="green")
#    ax4.plot(mat.sum(axis=0),color="royalblue")
#    
##    list_all_contact = list_all_contact+mat.sum(axis=0)
#    plt.xlabel("Position along the plasmid (bins 200 bp)")
#    plt.title("Contact of the plasmid with the host region")
#    
#    fig.tight_layout(pad=2.)
#    

    

# plot of sum:

plt.plot(list_all_contact)    
plt.xlabel("Position along the plasmid")
plt.ylabel("Contact with host hot spots of contact")













    