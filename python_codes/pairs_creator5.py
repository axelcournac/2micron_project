# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 18:15:46 2018
@author: axel KournaK 
To create all possible pairs from genomic positions.
"""
import numpy as np
import matplotlib
from pylab import *
import pandas as pd
import itertools

# Input: list of genomic positions :
df=pd.read_table('/home/axel/Bureau/YEAST/positions_peaks_condensins3.txt',header=None, delimiter=" ")
df=pd.read_table('/media/axel/RSG3/BACK_UP/data/yeasts_species_project/scerevisiae/bio_data/trna.dat2',header=None, delimiter="\t")
df=pd.read_table('/media/axel/RSG3/BACK_UP/data/yeasts_species_project/scerevisiae/bio_data/ori_precoces_Alvino.dat2',header=None, delimiter=" ")
df=pd.read_table('/media/axel/RSG3/BACK_UP/data/yeasts_species_project/scerevisiae/bio_data/rloops.txt2',header=None, delimiter="\t")
df=pd.read_table('/home/axel/Bureau/z_python_scripts/positions_HEG.txt',header=None, delimiter=" ")
df=pd.read_table('/home/axel/Bureau/z_python_scripts/positions_HEG_20pc.txt',header=None, delimiter=" ")
df=pd.read_table('/home/axel/Bureau/z_python_scripts/list_50pc_expr_coh_peaks.txt',header=None, delimiter=" ")
df=pd.read_table('list_20pc_expr_coh_peaks.txt',header=None, delimiter=" ")
df=pd.read_table('/home/axel/Bureau/z_python_scripts/list_50pc_expr_coh_peaks.txt',header=None, delimiter=" ")

df=pd.read_table('/media/axel/RSG41/diverse_yeast_data/CHip_seq_2018/test_summits.bed',header=None, delimiter="\t")

df=pd.read_table('/home/axel/Bureau/YEAST/Rloops_data/rloops_in_sup1.txt3',header=None, delimiter=" ")

df=pd.read_table('/media/axel/RSG4/diverse_yeast_data/ChIP_seq_quiescence/positions_Brn1_in_Q_All_Replicates_2000.txt',header=None, delimiter=" ")

df=pd.read_table('/home/axel/Bureau/YEAST/centro1.dat4',header=None, delimiter=" ")
df=pd.read_table('/home/axel/Bureau/YEAST/Rloops_data/rloops_sup_500bp.txt4',header=None, delimiter=" ")

df=pd.read_table('/home/axel/Bureau/YEAST/convergent_TSS_bins2kb.txt',header=None, delimiter=" ")
df=pd.read_table('/media/axel/RSG4/z_python_scripts/convergent_bins2kb.txt',header=None, delimiter=" ")

df=pd.read_table('/home/axel/Bureau/YEAST/groupe_common_polII_Log_reverse.txt',header=None, delimiter=" ")
df=pd.read_table('/home/axel/Bureau/YEAST/groupe_common_polII_Log_reverse.txt',header=None, delimiter=" ")


df=pd.read_table('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/positions_peaks_PolII_Log.txt',header=None, delimiter=" ")

df=pd.read_table('/home/axel/Bureau/YEAST/centro1.dat4',header=None, delimiter=" ")

df=pd.read_table('/home/axel/Bureau/YEAST/HSC_plasmids_in.txt2',header=None, delimiter="\t")
df=pd.read_table('/home/axel/Bureau/YEAST/GENES_SC288/long_genes_only_host_chrm.txt34',header=None, delimiter="\t")

df=pd.read_table('/home/axel/Bureau/figure_all_chrms_sup1/redone_norm3/th_08/HSC_plasmids_in_Micro-C_WT_log_SC288_genome.txt.sort.formated3.sorted',header=None, delimiter="\t")

# First implementation :                  
#   Output:  
#f_out = open("pairs_groupe_common_polII_Log_forward.txt","w+")
##BIN= 2000
#BIN=1
#
#list_all_chrms = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16") 
#
## intra pairs :
#for chr1 in list_all_chrms:
#    print(chr1)
#    bs = df.loc[(df[0] == chr1)]
#    bs = bs[1]
#    print(len(bs) )
#    combi_pos = list(itertools.combinations(bs , 2) )
#    for p in range(0, len(combi_pos) ):    #  to write pairs of positions of peaks and potential loops
#        f_out.write(chr +' '+ str( int(combi_pos[p][0]/BIN) ) + ' ' +str( int(combi_pos[p][1]/BIN))  +'\n')
#        
#f_out.close()

#------------------------------------------------------------------------------
# Intra and inter pairs :

f_out1 = open("pairs_intra_73HSC.txt","w+")
f_out2 = open("pairs_inter_73HSC.txt","w+")

len(df)
list_pairs = list(itertools.combinations( range(len(df)) ,2 ))
len(list_pairs)    

BIN=2000

for e in list_pairs :
    e1, e2 = e
    chr1, pos1 = df.iloc[e1]
    chr2, pos2 = df.iloc[e2]
    if chr1 == chr2 : 
        f_out1.write(chr1 + '\t'+ str( int(pos1) )+ '\t' + str( int(pos1+BIN))+ '\t' + 
                     chr1 + '\t'+ str( int(pos2)) + '\t' + str( int(pos2+BIN))+ '\t' + "1" +'\n')
    if chr1 != chr2:
        f_out2.write(chr1 + '\t'+ str( int(pos1) )+ '\t' + str( int(pos1+BIN))+ '\t' + 
                     chr2 + '\t'+ str( int(pos2)) + '\t' + str( int(pos2+BIN))+ '\t' + "1" +'\n')

f_out1.close()
f_out2.close()






# varisous tests:  

list_all_chrms = ("chr5") 
BIN=1

chr="chr5"
bs = df.loc[(df[0] == chr)]
bs = bs[1]
len(bs)
combi_pos = list(itertools.combinations(bs , 2) )
for p in range(0, len(combi_pos) ):    #  to write pairs of positions of peaks and potential loops
    plot(int(combi_pos[p][0]/BIN), int(combi_pos[p][1]/BIN),'o', color="yellow")
