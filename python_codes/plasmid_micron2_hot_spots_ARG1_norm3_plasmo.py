#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 17:08:57 2020
@author: axel
To extract plasmid p2 micron vs chrms from cool files. 
hot spots: to carefully extract hot spots of contact 
between plasmid and host chromosomes 
"""
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf
import sys
import os
from scipy.signal import find_peaks
from scipy import interpolate

import numpy as np
import cooler

sys.path.insert(0, os.path.abspath("/home/axel/Bureau/z_python_scripts_copy"))

# contact data:
cool_file = sys.argv[1]
plasmid_chosen = sys.argv[2]
name_bank = sys.argv[3]

# cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 /home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/sacCer3.chr_sizes.txt:2000 valid_idx_pcrfree.pairs valid_idx_pcrfree.pairs.cool
c = cooler.Cooler(cool_file)
cooler.balance_cooler(c, store=True, mad_max=5)   # Normalisation 

d=c.info
total_reads = d['sum']
# SC288 assembly:
BIN=int(sys.argv[6])
#list_chr =  ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
#             "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
#             "chr16","chrM","plasmid_p2-micron"]
list_chr1 = c.chromnames
#list_chr1 = ["NC_008285.1"]

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
centro_dist['pKan-STB-P'] = 0
centro_dist['pKan_DeltaSTB'] = 0
centro_dist['pKan_DeltaREP1'] = 0
centro_dist['KANMX'] = 0
centro_dist['Pf3D7_API_v3'] = 0 

print(file_size)
print(file_centro)
# function for nan sum
def nan_sum(row):
    total = np.nanmean(row)*len(row[~np.isnan(row)])
    return total

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
sizes_dist['pKan-STB-P'] = 12648
sizes_dist['pKan_DeltaSTB'] = 12313
sizes_dist['pKan_DeltaREP1'] = 11281
sizes_dist['KANMX'] = 1357
sizes_dist['Pf3D7_API_v3'] = 28798
#
genes = pd.read_csv('/home/axel/Bureau/YEAST/GENES_SC288/orf_genomic_all.txt2',sep=" ", header= None)
#genes = pd.read_csv('/home/axel/Bureau/walti/all_genes.txt3', 

long_genes = genes[ abs(genes[3]-genes[2]) > 5000 ]
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
limit_y =  0.00175
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

#list_chr1 =  ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
#             "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
#             "chr16","chr17","chr18","chr19","chr20","chr21","chr22",
#             "chrX"]

#list_chr1 =  ["chrX"]

f_out  = open("figures_per_chr.txt","w")
ci=0
for chr1 in list_chr1 :
    ci=ci+1
    long_genes_chr = long_genes[long_genes[0] ==chr1]
    ori_chr = ori[ori[0] ==chr1]
    ars_chr = ars[ars[0] ==chr1]
    acetyl_chr = acetyl[acetyl[0] ==chr1]

    # previous method
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
    print(chr1,reads_chr1,sum(coverage), sizes_dist[chr1], 
          reads_chr12/sizes_dist[chr1],reads_chr12/sum(coverage))
    
    f_out.write(chr1 +' '+ str( reads_chr1)+
                ' ' +str(sum(coverage) ) + 
                ' ' +str( sizes_dist[chr1]) +
                ' ' +str(reads_chr12/sizes_dist[chr1] ) +
                ' ' +str(reads_chr12/sum(coverage)) + '\n')
    
    coverage_plasmid=coverage_plasmid/coverage # normalisation 
    coverage_plasmid=coverage_plasmid/perc_plasmid # taking into accournt plasmid presence

    # computation of correlation:
    maxi= matscn.shape[0]*BIN
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

    col_chr = [chr1] * len(hot_spots[0])
    col_pos = list(hot_spots[0])
    col_pos = [int(e*BIN) for e in col_pos]
    df = pd.DataFrame({'chrom1': col_chr, 'start1': col_pos})
    df_total = pd.concat([df_total, df])


    # Multiplot:
    fig = plt.figure()
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(8, 12.8)
    fig.tight_layout(pad=5.0)
    gs = gridspec.GridSpec(4, 1,height_ratios=[7,1,1,1])

    ax1 = plt.subplot(gs[0])
    ax1.imshow(matscn**0.15, interpolation="none", cmap="afmhot_r")
    locs, labels = plt.xticks()            # Get locations and labels
    locs2 = [int(i*BIN) for i in locs]
    plt.xticks(locs,locs2)
    fig.tight_layout(pad=2.)
    plt.title("\n\n\n\n"+chr1+" "+str(total_reads)+" "
              +str(reads_chr1)+" "+str(reads_chr2)+" "+str(reads_chr12)+"  "
              +name_bank)

    ax2 = plt.subplot(gs[1], sharex=ax1)
    ax2.plot(coverage1)
    plt.title("Coverage")

#    coverage_plasmid = coverage_plasmid /(float(sizes_dist[chr2]))
    median_plasmid = np.nanmedian(coverage_plasmid)
    print(median_plasmid)
    ax3 = plt.subplot(gs[2], sharex=ax1)
    ax3.plot(coverage_plasmid_new, color="royalblue" )
    ax3.fill_between(range(len(coverage_plasmid_new)),coverage_plasmid_new, facecolor='royalblue')
    ax3.plot(float(centro_dist[chr1])/BIN,0.0,'o', color="red", markersize=4)
#    plt.xticks([])

    # plot with rectangles: 
    pos1=np.array(long_genes_chr[2])
    pos2=np.array(long_genes_chr[3])
    for i in range(len(pos1)):
        rectangle = plt.Rectangle( (min(pos1[i],pos2[i])/BIN,0.0),
                                  abs(pos2[i]-pos1[i])/BIN, 
                                  limit_y/4, 
                                  fc="grey")
        plt.gca().add_patch(rectangle)

    plt.axhline(y=median_plasmid,ls='--',color="Black")
#    plt.xticks([])
#    plt.ylim(0, limit_y)
#    plt.xlim(0, limit_x)
#    plt.plot(hot_spots[0],hot_spots[0]/hot_spots[0]*limit_y,'v', color="orange")
    h=hot_spots[1]
    left_ips = h['left_ips']
    right_ips = h['right_ips']
#    plt.plot(left_ips,left_ips*0.0,'|')
#    plt.plot(right_ips,right_ips*0.0,'|')
    plt.title(chr2+" Median="+str(round(median_plasmid*10**4,2))+"x 10^-4")

    ax4 = plt.subplot(gs[3], sharex=ax1)
    ax4.plot(coverage_plasmid_new, color="darkblue")

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

df_all.to_csv("contact_signal2_"+chr2+"_"+name_bank+".txt", sep='\t',
                index=False,na_rep='NA')

f_out.close()