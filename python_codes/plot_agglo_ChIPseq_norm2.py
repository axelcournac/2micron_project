#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 15:53:23 2021
@author: axel
"""
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


name_bank="Rpb1-1"
list_file=[""]
list_names=[""]


name_bank="tsuki2_paper"
list_file=["/media/axel/RSG4/H4_patch_mutants_tsuky/73_detected_loci/1D_ENRICHMENT_mean_Rpb3_in_Log_HHF2_Rep1_HSC_plasmid_in_log_recomputed_2000bp.txt",
           "/media/axel/RSG4/H4_patch_mutants_tsuky/73_detected_loci/1D_ENRICHMENT_mean_Rpb3_in_Log_HHF2_Rep2_HSC_plasmid_in_log_recomputed_2000bp.txt",
           "/media/axel/RSG4/H4_patch_mutants_tsuky/73_detected_loci/1D_ENRICHMENT_mean_Rpb3_in_Q_HHF2_Rep2_HSC_plasmid_in_log_recomputed_2000bp.txt",
           "/media/axel/RSG4/H4_patch_mutants_tsuky/73_detected_loci/1D_ENRICHMENT_mean_Rpb3_in_Q_5toA_Rep1_HSC_plasmid_in_log_recomputed_2000bp.txt",
           "/media/axel/RSG4/H4_patch_mutants_tsuky/73_detected_loci/1D_ENRICHMENT_mean_Rpb3_in_Q_R17R19A_Rep1_HSC_plasmid_in_log_recomputed_2000bp.txt"]
list_names=["Rpb3_in_Log_HHF2_Rep1","Rpb3_in_Log_HHF2_Rep2","Rpb3_in_Q_HHF2_Rep2",
            "Rpb3_in_Q_5toA_Rep1","Rpb3_in_Q_R17R19A_Rep1"]

name_bank="tsuki2_paper_3"
list_file=["/home/axel/Bureau/H4_patch_mutants_tsuky/73_detected_loci/1D_ENRICHMENT_mean_Rpb3_in_Log_HHF2_Rep1_HSC_plasmid_in_log_recomputed_2000bp.txt",
           "/home/axel/Bureau/H4_patch_mutants_tsuky/73_detected_loci/1D_ENRICHMENT_mean_Rpb3_in_Q_5toA_Rep1_HSC_plasmid_in_log_recomputed_2000bp.txt"
]
list_names=["Rpb3_in_log_HHF2_Rep1",
            "Rpb3_in_log_5toA_Rep1"]

name_bank="tsuki2_paper_3"
list_file=["/home/axel/Bureau/H4_patch_mutants_tsuky/73_detected_loci/1D_ENRICHMENT_mean_Rpb3_in_Q_HHF2_Rep2_HSC_plasmid_in_log_recomputed_2000bp.txt",
           "/home/axel/Bureau/H4_patch_mutants_tsuky/73_detected_loci/1D_ENRICHMENT_mean_Rpb3_in_Q_5toA_Rep2_HSC_plasmid_in_log_recomputed_2000bp.txt",
           "/home/axel/Bureau/H4_patch_mutants_tsuky/73_detected_loci/1D_ENRICHMENT_mean_Rpb3_in_Q_R17R19A_Rep2_HSC_plasmid_in_log_recomputed_2000bp.txt"]
list_names=["Rpb3_in_Q_HHF2_Rep2","Rpb3_in_Q_5toA_Rep1","Rpb3_in_Q_R17R19A_Rep1"]


# another coding of colours 
from colour import Color
pink = Color("pink")
royalblue = Color("royalblue")
colors = list(royalblue.range_to(Color("pink"), len(list_names) +1 ))

colors2=["red","pink","gray"] 
#-------------------------   
chr2="plasmid_p2-micron"   
BIN=2000
area=50000    # around the hot spots of contact of plasmid 
# Agglomerated plot: 
plt.title("Profile of ChIPseq "+chr2+"\nwith the "
          +" host regions of contact")
      
figure(figsize=(3, 2.7))
i=0
j=0
for file in list_file:
    i+=1
    color1=lighten_color('royalblue', i*0.2)
    color1= colors[i].hex_l
    f=np.loadtxt(file)
    moy = f[:,0]
    std = f[:,1]
    #moy = np.nanmedian(list_all_contact2,axis=0)
    #std = np.nanstd(list_all_contact2,axis=0)
    plt.plot(moy, linewidth=3., label=list_names[j], color=colors2[i-1])
    
    plt.fill_between(range(len(moy)),(moy-std/2.0),
                     (moy+std/2.0), color=colors2[i-1],alpha=.33)
    
    plt.vlines(x=area/BIN,ymin=0.0,ymax=max(moy)*1.5, 
               linestyles='--',
               color="black")
    
    plt.xlabel("Position around Hot spots of contact with 2u plasmid")
    plt.ylabel("ChIPseq / Input")
    
    tick_locs=(0,int(area/BIN),int(2*area/BIN))
    tick_lbls=('-'+str(area/BIN)+' kb','Set','+'+str(area/1000)+' kb')
    plt.xticks(tick_locs, tick_lbls,rotation=45)
    
    plt.subplots_adjust(bottom=0.2)
    plt.subplots_adjust(left=0.2)
    plt.ylim(0.0,1.8)
#    plt.xlim(0.,area*2/BIN)
    plt.xlim(10,area*2/BIN-10)
    j+=1

plt.legend(loc = "lower left")
#plt.legend(loc = "upper right")

plt.savefig("agglomerated_signal_on_HSC_"+"_"+name_bank+"_"+str(BIN/1000)+'.pdf')
plt.close("all")

moy=moy.reshape(int(area/BIN*2+1),1)
std=std.reshape(int(area/BIN*2+1),1)
c2= np.concatenate((moy, std), axis=1)

np.savetxt(name_bank+"_agglomerated_signal_on_HSC_"+
                "_"+name_bank+"_"+str(BIN/1000)+"kb_norm2"+'.txt', c2)

