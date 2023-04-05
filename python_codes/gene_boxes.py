#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 12:31:01 2023
@author: axel
to visualise gene and one signal 
2: we add color gene in fction of its size 
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random 

import matplotlib.colors as mc
import colorsys

genes = pd.read_csv('/home/axel/Bureau/YEAST/GENES_SC288/orf_genomic_all.txt2', 
                    sep=" ", header= None)
long_genes = genes[ abs(genes[3]-genes[2]) > 7000 ]
len(long_genes)

np.min(abs(genes[3]-genes[2]))
np.max(abs(genes[3]-genes[2]))

# 1)
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

# 2) Another coding of colours 
from colour import Color
white = Color("cyan")
royalblue = Color("blue")
colors = list(white.range_to(Color("blue"), 12 +1 ))

# colorbar 
for i in range(6):
    print(i)
    j=i
    if j>=5:
        j=5
#     1)
    color1=lighten_color('darkblue', j/5.0)
#     2)
#    color1= colors[i].hex_l
    
    plt.plot(i,0,marker="s",color=color1, markersize=50.0)

      
window = 20000
file_HSC='/home/axel/Bureau/figure_all_chrms_sup1/redone_norm3/th_08/HSC_plasmids_in_Micro-C_WT_log_SC288_genome.txt.sort'
df= pd.read_table(file_HSC,header=0,sep="\t")

chr1="chr11"

bin_size=1000 
for ir in range(len(df)):                 
    chr1=df['chrom1'][ir]
    pos=int(df['pos'][ir])
    
    pos1=pos-window
    pos2=pos+window
    
    print(chr1+"\t"+str(pos-5)+"\t"+str(pos+5))
    
    genes_chr = genes[genes[0] ==chr1]
    
    # plot of genes:
    genes_chr=np.array(genes_chr)
    for j in range(len(genes_chr)) :
        y_rand=random.uniform(0, 1)
        y_rand=0
        p1=int(genes_chr[j][2])
        delta_p=int( (genes_chr[j][3]-genes_chr[j][2]))
        
        size_a=np.abs(delta_p)
        size_i=np.int(size_a/bin_size)
        if size_i>=5:
            size_i=5
        if size_i==0:
            size_i=1
        # 2 ways of coloring 
        color1=lighten_color('darkblue', size_i/5.0)
#        color1= colors[size_i].hex_l
        
        plt.arrow(p1, y_rand, delta_p, 0,
                  width=5,length_includes_head=True,head_length=200,
                  head_width=10, color=color1) 
    #        plt.text(p1+delta_p/2.0,y_rand,genes_chr[j][1]) 
    plt.tight_layout()    
       
    fig = plt.gcf()
    fig.set_size_inches(18.5, 0.5)
    plt.tight_layout()
    plt.xlim(pos1,pos2)
    plt.axis('off')
    
    fig.savefig('window7_gene_'+str(ir)+"_"+chr1+'_'+str(pos)+'.pdf',dpi=400)
    plt.close("all")  
    #plt.xlim(530000,555000)    
    #plt.ylim(-2,2)
    
  # pdfjam $(ls -v  window7_gene_*.pdf) --nup 1x73 --landscape --outfile Page12_win_genes7.pdf 
    
    
    