#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 17:16:09 2024
@author: axel
To have the expression of yeast genes and 2µ genes from tables 
"""

import pandas as pd
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.ticker as mtick

# input aligned with bowtie2 and counts made with htseq-count software 

# file_rna='/home/axel/Bureau/transcription_level_plasmid_genes/gene_count_ERR12074422_with_2u.txt'
# file_rna='/home/axel/Bureau/transcription_level_plasmid_genes/gene_count_FG0222_with_2u.txt'
file_rna='/home/axel/Bureau/transcription_level_plasmid_genes/gene_count_SRR13736618_IP_with_2u.txt'
rna = pd.read_csv(file_rna, sep="\t", header= None, skipfooter=5, engine='python')


# file_shot='/home/axel/Bureau/transcription_level_plasmid_genes/gene_count_ERR12074416_with_2u.txt'
# file_shot='/home/axel/Bureau/transcription_level_plasmid_genes/gene_count_FG0223_with_2u.txt'
file_shot='/home/axel/Bureau/transcription_level_plasmid_genes/gene_count_SRR13736616_input_with_2u.txt'

shot = pd.read_csv(file_shot, sep="\t", header= None, skipfooter=5, engine='python')

# Normaliosation step: we divide by total sum of read count
rna[1]=np.divide(rna[1],np.sum(rna[1]), where=True, dtype='float64')
shot[1]=np.divide(shot[1],np.sum(shot[1]), where=True, dtype='float64')

np.set_printoptions(precision=10, floatmode='maxprec')

#------------------------------------------------------------------------------
# simple plot : 
plt.plot(shot[1],rna[1]/shot[1],'o')

rna_FLPa= rna[rna[0]=="FLPa"]
shot_FLPa= shot[shot[0]=="FLPa"]
plt.plot(shot_FLPa[1], rna_FLPa[1]/shot_FLPa[1],'o')
plt.text(shot_FLPa[1], rna_FLPa[1]/shot_FLPa[1], 'FLPa')

rna_FLPb= rna[rna[0]=="FLPb"]
shot_FLPb= shot[shot[0]=="FLPb"]
plt.plot(shot_FLPb[1], rna_FLPb[1]/shot_FLPb[1],'o')
plt.text(shot_FLPb[1], rna_FLPb[1]/shot_FLPb[1], 'FLPb')

rna_REP1= rna[rna[0]=="REP1"]
shot_REP1= shot[shot[0]=="REP1"]
plt.plot(shot_REP1[1], rna_REP1[1]/shot_REP1[1],'o')
plt.text(shot_REP1[1], rna_REP1[1]/shot_REP1[1],'REP1')

rna_REP2= rna[rna[0]=="REP2"]
shot_REP2= shot[shot[0]=="REP2"]
plt.plot(shot_REP2[1], rna_REP2[1]/shot_REP2[1],'o')
plt.text(shot_REP2[1], rna_REP2[1]/shot_REP2[1],'REP2')

rna_RAF= rna[rna[0]=="Raf"]
shot_RAF= shot[shot[0]=="Raf"]
plt.plot(shot_RAF[1], rna_RAF[1]/shot_RAF[1],'o') 
plt.text(shot_RAF[1], rna_RAF[1]/shot_RAF[1],'Raf') 

plt.xlabel("Number of reads from DNA signal for each gene")

plt.ylabel("Transcription level = Number of reads from RNA-seq signal / Number of reads from DNA signal for each gene")


#------------------------------------------------------------------------------
# With histogram included:
x1=shot[1]
y1=rna[1]/shot[1]

x1=np.log10(x1)
y1=np.log10(y1)

x1 = x1[np.isfinite(y1)]
y1 = y1[np.isfinite(y1)]


fig = plt.figure(figsize=(12,5))
gs = gridspec.GridSpec(6, 6)
ax_main = plt.subplot(gs[0:6, :5])
#ax_xDist = plt.subplot(gs[0, :2],sharex=ax_main)
ax_yDist = plt.subplot(gs[0:6, 5],sharey=ax_main)
plt.setp(ax_yDist.get_yticklabels(), visible=False)  
  
ax_main.scatter(x1,y1,marker='o', alpha=.5)
ax_main.set_title(file_rna)
ax_main.set(xlabel=r"$\bf{DNA}$"+" "+r"$\bf{content}$"+ " = Number of reads from DNA signal for each gene (input, log10)",
            ylabel=r"$\bf{Transcription}$"+" "+r"$\bf{Level}$"+ "\n = Number of reads from ChIP-seq Brn1 signal \n divided by Number of reads from DNA signal for each gene")


rna_FLPa= rna[rna[0]=="FLPa"]
shot_FLPa= shot[shot[0]=="FLPa"]
ax_main.plot(np.log10(shot_FLPa[1]), np.log10(rna_FLPa[1]/shot_FLPa[1]),'o')
ax_main.text(np.log10(shot_FLPa[1]), np.log10(rna_FLPa[1]/shot_FLPa[1]), 'FLPa')

rna_FLPb= rna[rna[0]=="FLPb"]
shot_FLPb= shot[shot[0]=="FLPb"]
ax_main.plot(np.log10(shot_FLPb[1]), np.log10(rna_FLPb[1]/shot_FLPb[1]),'o')
ax_main.text(np.log10(shot_FLPb[1]), np.log10(rna_FLPb[1]/shot_FLPb[1]), 'FLPb')

rna_REP1= rna[rna[0]=="REP1"]
shot_REP1= shot[shot[0]=="REP1"]
ax_main.plot(np.log10(shot_REP1[1]), np.log10(rna_REP1[1]/shot_REP1[1]),'o')
ax_main.text(np.log10(shot_REP1[1]), np.log10(rna_REP1[1]/shot_REP1[1]),'REP1')

rna_REP2= rna[rna[0]=="REP2"]
shot_REP2= shot[shot[0]=="REP2"]
ax_main.plot(np.log10(shot_REP2[1]), np.log10(rna_REP2[1]/shot_REP2[1]),'o')
ax_main.text(np.log10(shot_REP2[1]), np.log10(rna_REP2[1]/shot_REP2[1]),'REP2')

rna_RAF= rna[rna[0]=="Raf"]
shot_RAF= shot[shot[0]=="Raf"]
ax_main.plot(np.log10(shot_RAF[1]), np.log10(rna_RAF[1]/shot_RAF[1]),'o') 
ax_main.text(np.log10(shot_RAF[1]), np.log10(rna_RAF[1]/shot_RAF[1]),'Raf') 

ax_main.set_ylim(min(y1)*1.2,max(y1)*1.1)
ax_yDist.hist(y1,bins=50,orientation='horizontal',align='mid', color="grey")
ax_yDist.axhline(y = np.median(y1), color = 'r')
ax_yDist.set(xlabel='Nber of genes')
ax_yDist.set_ylabel('Transcription Level (log10)', loc='center')
ax_yDist.yaxis.set_label_position("right")
ax_yDist.yaxis.tick_right()

# #ax_yDist.ylim(min(),max())
# plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
# 




