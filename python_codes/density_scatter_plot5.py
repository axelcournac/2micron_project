# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 11:30:00 2016
@author: axel
"""
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec 
from scipy.stats import gaussian_kde
import pandas as pd
import re
import matplotlib.ticker as mtick
%matplotlib qt

#  simple scatter plot: 
#file="/home/axel/Bureau/2micron_plasmid_PROJECT/analysis_chip_exo_scatter/scores_all_proteins_HSC_plasmid_ratio_mean_2000bp.txt.sorted"
file="/home/axel/Bureau/2micron_plasmid_PROJECT/analysis_chip_exo_scatter/nov2021/scores_all_proteins_73HSC_plasmid_ratio_mean_2000bp.txt"
scatter = pd.read_table(file,header=None, delimiter=" ") 
print(len(scatter))
scatter= scatter.loc[ (scatter[3] > 500000) ]
print(len(scatter))

# dictionary :
names = pd.read_table("/home/axel/Bureau/2micron_plasmid_PROJECT/analysis_chip_exo_scatter/SraRunTable_ChIP_exo.txt",
                       header=None, delimiter=" ") 
zipbObj = zip(names[0], names[1])
# Create a dictionary from zip object
names_dict = dict(zipbObj)
z=  names[0]
z2 = [i[0:11] for i in z]
z3 = [ names_dict[i] for i in z2]
z4= [re.findall(r"^\d+_(\w+)_\w+_?-?ChIP", i) for i in z3]
zipbObj = zip(names[0], z4)
names2_dict = dict(zipbObj)

# dictionary for UMAP categories:
names2 = pd.read_table("/home/axel/Bureau/2micron_plasmid_PROJECT/analysis_chip_exo_scatter/UMAP_results_YEP/UMAP_categories.txt",
                       header=None, delimiter="\t") 
zipbObj = zip(names2[0], names2[9])
names2_dict = dict(zipbObj)

# With colors depending on categories from UMAP analysis:
cat = set(names2_dict.values())
colours = ("black","olive","orange","yellow","cyan","gold","purple","gray","wheat","seagreen",
     "dodgerblue","peru","aqua","plum","magenta","lime","brown","indigo","lightgreen")

zipbObj = zip(cat, colours)
color_dict = dict(zipbObj)

x=np.array(scatter[3])
y=np.array(scatter[2])
z=np.array(scatter[1])

for s in range(len(scatter)) :
    if z[s] in names2_dict.keys() :
        for n in cat:
            if names2_dict[z[s]] == n:
                plt.plot(x[s], y[s], 'o', color=color_dict[n])

# to put legend 
for n in cat:
    plt.plot(1000,0.4,'o', color=color_dict[n],label=n)
plt.plot(1000,0.4, 'o', alpha=0.5, label="Without category", color="blue")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.close('all')

# Plot sorted by category:
ic=0
th_sig=1.25  # Threshold to plot the name of the protein
for n in cat:
    print(n)
    for s in range(len(scatter)) :
        if z[s] in names2_dict.keys() :
            if names2_dict[z[s]] == n:
                ic=ic+1
                plt.plot(ic, y[s], 'o', color=color_dict[n])
                if y[s] > th_sig:
                    plt.text(ic, y[s], z[s])

for s in range(len(scatter)) :
    if z[s] not in names2_dict.keys() :
        ic=ic+1
        plt.plot(ic, y[s], 'o', color="blue", alpha=0.5)
        if y[s] > th_sig:
            plt.text(ic, y[s], z[s])

plt.xlabel("ChIP-exo library sorted by UMAP categories")
plt.ylabel("Average signal of ChIP-exo \n around Hot spots of contact of plasmid")

# to put legend 
for n in cat:
    plt.plot(1000,0.4,'o', color=color_dict[n],label=n)
plt.plot(1000,0.4, 'o', alpha=0.5, label="Without category", color="blue")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))