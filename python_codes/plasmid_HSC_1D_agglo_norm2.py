#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 17:08:57 2020
@author: axel
To extract plasmid p2 micron vs chrms from cool files.
To generate plasmid map and hot spot of contact on host genome
ALL: all hot spots of contact between plasmid and host 
To have the agglomerated signal around HSC made by a plasmid. 
"""
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf
import random
import sys
import os
sys.path.insert(0, os.path.abspath("/home/axel/Bureau/z_python_scripts_copy"))
import numpy as np
import cooler
from Bio import SeqIO

number_of_colors = 18
color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]

cool_file = sys.argv[1]
plasmid_chosen = sys.argv[2]
name_bank = sys.argv[3]
file_set_positions= sys.argv[4]

#----------
if not os.path.exists(name_bank+"_files") :
    os.makedirs(name_bank+"_files")

c = cooler.Cooler(cool_file)
cooler.balance_cooler(c, store=True, mad_max=10)   # Normalisation 

list_chr1 = c.chromnames

#file_set_positions= "/home/axel/Bureau/figure_all_chrms_sup1/redone/HSC_plasmids_in_Micro-C_WT_log_SC288_genome.txt.sort3"
# Regions where to look the contacts: 
#file_set_positions='/home/axel/Bureau/YEAST/HSC_plasmids_in.txt'
#pos_set = pd.read_table('/home/axel/Bureau/YEAST/strong_loops_cohesin_quant_cohe_classic_10-50kb_49_mitotic.txt')
pos_set = pd.read_table(file_set_positions)
#pos_set = pd.read_table('/home/axel/Bureau/YEAST/centro1_SC288.dat4',sep=" ")        
#pos_set = pd.read_table('/home/axel/Bureau/YEAST/HSC_plasmids_in_200pb.txt')
              
input_file_type = 1
name_pos_set= name_bank

if plasmid_chosen == "None":
    chr2 = "plasmid_p2-micron"   
else :
    chr2 = plasmid_chosen 

total_number_reads=c.info['sum']
reads_plasmide=0
for chr3 in list_chr1:
            m1 = c.matrix(balance=False).fetch(chr2, chr3)
            if chr2==chr3:
                reads_plasmide += m1.sum()/2+np.trace(m1)/2
            else:
                reads_plasmide += m1.sum()

perc_plasmid=reads_plasmide/total_number_reads
print("Percentage plasmide")
print(perc_plasmid)

#
BIN=2000
area=40000    # around the hot spots of contact of plasmid 
len_area_contact = int(area/BIN) *2 +1

if input_file_type == 1 :  # only one position to do the agglomeration
   d = {'chrom1': pos_set.chrom1, 'start1': pos_set.start1}
   p1 = pd.DataFrame(data=d)
   p12 = p1
   p12 = p12.drop_duplicates(subset=None, keep='first', inplace=False)

if input_file_type == 2 :  # some pairs positions to do the agglo
    d = {'chrom1': pos_set.chrom1, 'start1': pos_set.start1}
    p1 = pd.DataFrame(data=d)
    d = {'chrom1': pos_set.chrom2, 'start1': pos_set.start2}
    p2 = pd.DataFrame(data=d)
    
    p12  =  pd.concat([p1, p2])
    p12 = p12.drop_duplicates(subset=None, keep='first', inplace=False)
 
print(len(p12))

# SC288 assembly:
if chr2 == "plasmid_p2-micron":
    plasmid_genes = pd.read_csv('/home/axel/Bureau/YEAST/GENES_SC288/genes_plasmid_J01347.1.txt', 
                    sep=" ", header= None)
if chr2 ==  "pRS413" :
    plasmid_genes = pd.read_csv('/home/axel/Bureau/YEAST/genes_pRS413.txt3',
                            sep="\t", header= None)
    
centro = pd.read_table("/home/axel/Bureau/YEAST/centro1.dat4",
                       header=None, delimiter=" ") 
# Create a zip object from two lists
zipbObj = zip(centro[0], centro[1])
# Create a dictionary from zip object
centro_dist = dict(zipbObj)
centro_dist['chrM'] = 0
centro_dist['plasmid_p2-micron'] = 0
centro_dist['pJL638'] = 0
centro_dist['pYADE4'] = 0
centro_dist['pKan-STB-P'] = 0
centro_dist['pKan_DeltaSTB'] = 0
centro_dist['pKan_DeltaREP1'] = 0
centro_dist['pST1760'] = 0
centro_dist['pYES2'] = 0
centro_dist['pRS416'] = 0
centro_dist['KANMX'] = 0

sizes = pd.read_table("/home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/sacCer3.chr_sizes.txt",
                      header=None, delimiter="\t") 
zipbObj = zip(sizes[0], sizes[1])
sizes_dist = dict(zipbObj)
sizes_dist['plasmid_p2-micron']=6300
sizes_dist['pRS413'] = 4970  # the one of Fabien 
sizes_dist['pSH47'] = 6979
sizes_dist['pJL638'] = 4970
sizes_dist['pYADE4'] = 6037
sizes_dist['pKan-STB-P'] = 12648
sizes_dist['pKan_DeltaSTB'] = 12313
sizes_dist['pKan_DeltaREP1'] = 11281
sizes_dist['pRDN'] = 35913
sizes_dist['pLM006'] = 8377
sizes_dist['pST1760'] = 11000
sizes_dist['pYES2'] = 6363
sizes_dist['pRS416'] = 4898
sizes_dist['KANMX'] = 1357

genes = pd.read_csv('/home/axel/Bureau/YEAST/GENES_SC288/orf_genomic_all.txt2', 
                    sep=" ", header= None)
long_genes = genes[ abs(genes[3]-genes[2]) > 7000 ]
len(long_genes)

def nan_sum(row):
    total = np.nanmean(row)*len(row[~np.isnan(row)])
    return total
#---------
sizes_chr=c.chromsizes
repo_fasta="/home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/"
#repo_fasta="/home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron_without_ade2/"
for i in range(len(p12)) :
    print(i)
    chr1 = p12.chrom1[i]
    fasta_file=repo_fasta+chr1+".fa"
        
    record = SeqIO.read(open(fasta_file), "fasta")
    
    pos1 = p12.start1[i]-area + 1
    pos2 = p12.start1[i]+area + 1
    
    if pos1<=0:
        pos1=0
    if pos2>=len(record):
        pos2=len(record)
        
    if pos1<=sizes_chr[chr1] and pos2<=sizes_chr[chr1]:   
        chr1_focus = (chr1,pos1,pos2) # region of focus  
               
        # for the normalisation 2 
        matscn = c.matrix(balance=True).fetch(chr1, chr1)
        coverage1 = np.apply_along_axis(nan_sum, 0, matscn)
        mat = c.matrix(balance=True).fetch(chr1, chr2)
        coverage_plasmid1 = np.apply_along_axis(nan_sum, 1, mat)
    
        # to work with the raw number of reads 
        reads_chr1=0
        coverage=0
        for chr3 in list_chr1:
            m1 = c.matrix(balance=False).fetch(chr1, chr3)
            reads_chr1 += m1.sum()
            coverage_temp=np.apply_along_axis(nan_sum, 1, m1)
            coverage_temp[np.isnan(coverage1)] = np.nan
            coverage += coverage_temp
    
        m2 = c.matrix(balance=False).fetch(chr2, chr2)
        reads_chr2= m2.sum()
        contact_of_plasmid= np.apply_along_axis(nan_sum, 1, m2)
        m12 = c.matrix(balance=False).fetch(chr1, chr2)
        reads_chr12= m12.sum()
        coverage_plasmid = np.apply_along_axis(nan_sum, 1, m12)
        coverage_plasmid [ np.isnan(coverage_plasmid1)] = np.nan
    
        contact_with_plasmid=coverage_plasmid/coverage # norm2 here
        contact_with_plasmid=contact_with_plasmid/perc_plasmid # taking into accournt plasmid presence
        
        p1_binned = int(pos1/BIN)
        p2_binned = int(pos2/BIN)
        contact_with_plasmid = contact_with_plasmid[range(p1_binned,p2_binned+1)]
        
        # computation of correlation:
        maxi= matscn.shape[0]*BIN
        
        if pos1<=0:
            l1=len_area_contact-len(contact_with_plasmid)
            vect_nan = np.empty(l1)
            vect_nan[:] = np.NaN
            contact_with_plasmid  = np.concatenate((vect_nan, contact_with_plasmid))
            
        if pos2>=len(record):
            l1=len_area_contact-len(contact_with_plasmid)
            vect_nan = np.empty(l1)
            vect_nan[:] = np.NaN
            contact_with_plasmid  = np.concatenate((contact_with_plasmid, vect_nan))
            
    #    coverage[coverage==0] = np.nan
        if i==0:    
            list_all_contact = contact_of_plasmid
            list_all_contact2 = contact_with_plasmid 
        else : 
            list_all_contact = np.vstack((list_all_contact, 
                                                contact_of_plasmid))
            list_all_contact2 = np.vstack((list_all_contact2, 
                                                contact_with_plasmid))
    #   previous method:            
        matscn = c.matrix(balance=True).fetch(chr1_focus, chr2) # ! we focus on the region of contact
        matscn[np.isnan(matscn)] = 0 
        coverage = matscn.sum(axis=0)
        
        contact_with_plasmid = matscn.sum(axis=1)
        contact_of_plasmid = matscn.sum(axis=0)
    
        mat = c.matrix(balance=True).fetch(chr2, chr2)
        mat[np.isnan(mat)] = 0 
        plasmid_intra = mat.sum(axis=1)
    #    plasmid_intra[plasmid_intra==0] = np.nan
        
        #---------------------
        # Multiplot: 
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(8, 11.8)
        fig.tight_layout(pad=2.5)
        gs = gridspec.GridSpec(4, 2,height_ratios=[9,1,1,1],
                               width_ratios=[9,2])
    
        ax1 = plt.subplot(gs[0])
        ax1.imshow(matscn**0.5,interpolation="none", cmap="afmhot_r", 
                 aspect="auto")
        plt.ylabel(chr1+' '+str(int(p12.start1[i])))
        plt.xlabel(chr2)
        plt.title(chr1)
        
        ax11 = plt.subplot(gs[1],sharey=ax1)  # vertical
        v1 = range(matscn.shape[0])
        v2 = matscn.sum(axis=1)
    #    plt.plot(0.0,area/BIN,'^',color="orange")
        plt.hlines(y=area/BIN,xmin=0.0, xmax=max(v2), linestyles='--', label='')
        plt.title("Contacts with \n"+chr2)
        ax11.plot(v2,v1)
        plt.xticks(rotation=70)
        
        ax2 = plt.subplot(gs[2], sharex=ax1)
        plt.title("Proportion of the signal in INTRA of the plasmid (in %)") 
        plasmid_intra = plasmid_intra /(sizes_dist[chr2])
        median_plasmid = np.nanmedian(plasmid_intra) 
        ax2.plot(plasmid_intra)
        print(median_plasmid)
        
        ax3 = plt.subplot(gs[4], sharex=ax1)
        plt.ylim(-0.25, 0.25)
        plt.title("Plasmid Genes")
    #    plt.title(chr2+" Median="+str(round(median_plasmid*10**8,2))+"x 10^8")
    
        ax4 = plt.subplot(gs[6], sharex=ax1)
    #    ax4.plot(coverage_chrM, color="green")
        ax4.plot(matscn.sum(axis=0),color="royalblue")
        
        plt.xlabel("Position along the plasmid (bins 200 bp)")
        plt.title("Contact of the plasmid with the host region")
        
        fig.tight_layout(pad=2.)
        
        plt.savefig(name_bank+"_files/"+'_'+chr1+'_'+str(pos1)+
                    "_"+name_bank+"_"+str(BIN/1000)+"kb_scn"+'.pdf')
plt.close("all")    

#  Plottings:  
plt.imshow(list_all_contact2, interpolation="none")
plt.savefig(name_bank+"_files/all_HSC_"+
                "_"+name_bank+"_"+str(BIN/1000)+"kb_scn"+'.pdf')
plt.close("all")    
   
# Agglomerated plot: 
plt.title("Profile of contacts of plasmid "+chr2+"\nwith the "
          +str(len(p12))+" HSC, perc. of plasmid "+str(perc_plasmid))

moy = np.nanmedian(list_all_contact2,axis=0)
std = np.nanstd(list_all_contact2,axis=0)

plt.plot(moy, linewidth=3., color="royalblue", 
         label=chr2+" with "+name_pos_set)

plt.fill_between(range(len(moy)),(moy-std/2.0),
                 (moy+std/2.0),alpha=.33, 
                 color = 'royalblue')

plt.vlines(x=area/BIN,ymin=0.0,ymax=max(moy), 
           linestyles='--',
           color="black")

plt.xlabel("Position around "+name_pos_set)
plt.ylabel("Contact with plasmid "+chr2)

tick_locs=(0,int(area/BIN),int(2*area/BIN))
tick_lbls=('-'+str(area/1000)+' kb','Set','+'+str(area/1000)+' kb')
plt.xticks(tick_locs, tick_lbls,fontsize=15,rotation=45)

plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
#plt.ylim(0.0,0.0011)
plt.xlim(0.,40.)
plt.savefig(name_bank+"_files/agglomerated_signal_on_HSC_"+
                "_"+name_bank+"_"+str(BIN/1000)+"kb_scn"+'.pdf')
plt.close("all")

moy=moy.reshape(41,1)
std=std.reshape(41,1)
c2= np.concatenate((moy, std), axis=1)

np.savetxt(name_bank+"_files/agglomerated_signal_on_HSC_"+
                "_"+name_bank+"_"+str(BIN/1000)+"kb_norm2"+'.txt', c2)
