#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 15:53:23 2021
@author: axel
_n: generalised for n conditions 
"""
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
import sys 

n= int(sys.argv[1])

list_file=[]
list_names=[]
name_bank="_"

for i in range(n):
    list_file.append(sys.argv[i+2])
    list_names.append(sys.argv[n+i+2])
    name_bank=name_bank+"_"+sys.argv[n+i+2]

name1= sys.argv[3]
name2= sys.argv[4]

# if we need some nuances of color (for cell cycle)
def lighten_color(color, amount=0.5) :
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

#name_bank=name1+"_"+name2
#list_file=[file1, file2]
#list_names=[name1, name2]
       
# another coding of colours 
from colour import Color
pink = Color("pink")
royalblue = Color("royalblue")
colors = list(royalblue.range_to(Color("pink"), len(list_names) +1 ))

#-------------------------   
chr2="plasmid_p2-micron"   
BIN=2000
area=40000    # around the hot spots of contact of plasmid 
# Agglomerated plot: 
plt.title("Profile of contacts of plasmid "+chr2+"\nwith the "
          +" host regions of contact")
      
figure(figsize=(3, 2.7))
i=0
j=0
for f in list_file:
    i+=1
    color1=lighten_color('royalblue', i*0.2)
    color1= colors[i].hex_l
    f=np.loadtxt(f)
    moy = f[:,0]
    std = f[:,1]

    plt.plot(moy, linewidth=3., label=list_names[j])
    
    plt.fill_between(range(len(moy)),(moy-std/2.0),
                     (moy+std/2.0),alpha=.33)
    
#    plt.plot(moy, linewidth=3., label=list_names[j], color=color1)
#    
#    plt.fill_between(range(len(moy)),(moy-std/2.0),
#                     (moy+std/2.0),alpha=.33, color=color1)
    
    plt.vlines(x=area/BIN,ymin=0.0,ymax=max(moy), 
               linestyles='--',
               color="white", alpha=0.0)
    
#    plt.xlabel("Position around Hot spots of contact with 2u plasmid")
#    plt.ylabel("Contact with plasmid "+chr2)
    
    tick_locs=(0,int(area/BIN),int(2*area/BIN))
    tick_lbls=('-'+str(area/1000)+' kb','Set','+'+str(area/1000)+' kb')
    plt.xticks(tick_locs, tick_lbls,rotation=0)
    
    plt.locator_params(axis='y', nbins=3)
    
    plt.subplots_adjust(bottom=0.2)
    plt.subplots_adjust(left=0.2)
#    plt.ylim(0.0,0.8)
    plt.xlim(0.,area*2/BIN)
    j+=1

plt.legend(loc = "lower left")
#plt.legend(loc = "upper right")
plt.ylim(0,1.15)

plt.savefig("agglomerated_signal_on_HSC2_"+"_"+name_bank+"_"+str(BIN/1000)+'.pdf')
plt.close("all")

moy=moy.reshape(41,1)
std=std.reshape(41,1)
c2= np.concatenate((moy, std), axis=1)

np.savetxt(name_bank+"_agglomerated_signal_on_HSC_"+
                "_"+name_bank+"_"+str(BIN/1000)+"kb_norm2"+'.txt', c2)

