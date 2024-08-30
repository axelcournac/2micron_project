#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 13:33:39 2023
@author: axel
To compute overlab between contacted regions like peaks 
and biological objects like TY (Transposons of Yeast), long genes etc 
or sub-telomeric regions 
"""

import pandas as pd
import pyranges as pr
import numpy as np
import matplotlib.pyplot as plt

# conversion of the peaks position from Fab
# dicty
df=  {'chr01': [20, 33, 103, 175, 251, 345, 408, 442, 559, 595, 748, 805, 845, 936, 960, 1035, 1047, 1080, 1120, 1186, 1228, 1254, 1261, 1269, 1308, 1366, 1460, 1484, 1504, 1558, 1615, 1651, 1716, 1750, 1809, 1827, 1935, 1964, 2000, 2018, 2168, 2236, 2293, 2452], 
'chr02': [25, 51, 104, 133, 171, 294, 373, 428, 582, 621, 659, 727, 779, 936, 974, 1036, 1071, 1097, 2016, 2078, 2135, 2345, 2379, 2519, 2567, 2600, 2644, 2737, 2753, 2822, 2908, 2971, 2995, 3068, 3125, 3323, 3360, 3390, 3418, 3500, 3549, 3834, 3868, 4056, 4101, 4110, 4122, 4182, 4203, 4220], 
'chr03': [11, 61, 113, 205, 228, 377, 422, 554, 620, 655, 847, 976, 1014, 1080, 1145, 1253, 1315, 1330, 1423, 1514, 1637, 1708, 1748, 1800, 1836, 1874, 1905, 1946, 1965, 1994, 2010, 2043, 2138, 2202, 2218, 2369, 2558, 2662, 2673, 2707, 2767, 2863, 2967], 
'chr04': [163, 200, 265, 405, 493, 517, 663, 705, 840, 884, 905, 1054, 1280, 1467, 1512, 1658, 1698, 1763, 1828, 1861, 1924, 2071, 2091, 2218, 2301, 2420, 2503, 2553, 2664], 
'chr05': [10, 19, 55, 138, 189, 219, 253, 405, 531, 689, 771, 843, 1033, 1092, 1184, 1243, 1268, 1306, 1373, 1435, 1472, 1548, 1605, 1645, 1657, 1686, 1729, 1767, 1794, 1842, 1930, 2004, 2047, 2070, 2109, 2199, 2230, 2239, 2265, 2315, 2356, 2393, 2472, 2493, 2523], 
'chr06': [10, 81, 110, 136, 178, 199, 224, 272, 288, 370, 412, 511, 587, 732, 804, 863, 874, 899, 921, 989, 1021, 1149, 1171, 1193, 1263, 1310, 1367, 1416, 1434, 1616, 1677, 1720, 1771]}   
    

# Fermentati
df= {'chr01': [23, 131, 171, 368],
 'chr02': [11, 276, 388],
 'chr03': [40, 136, 190, 241, 373, 440, 477, 522],
 'chr04': [174, 254, 274, 309, 412, 460, 500],
 'chr05': [187, 220, 323, 525, 599],
 'chr06': [26, 226, 327, 425, 501, 563],
 'chr07': [41, 96, 147, 242, 355, 611, 696, 744],
 'chr08': [27, 240, 366, 383, 421, 451, 520, 578, 644, 711, 779]}

# waltii 
df={'chr01': [184, 244, 332, 401, 535],
 'chr02': [23, 104, 172, 220, 235, 402, 509],
 'chr03': [96, 145, 238, 388, 419, 504],
 'chr04': [8, 40, 121, 159, 264, 299, 390, 411, 477, 586],
 'chr05': [61, 115, 168, 309, 422, 461, 528, 589, 640],
 'chr06': [50, 151, 202, 229, 517, 624, 660],
 'chr07': [100, 153, 222, 268, 425, 582, 632],
 'chr08': [10, 260, 379, 483]}

for k in df.keys():
    # print(k)
    for pos in df[k]:
        print(k+"\t"+str(pos*2000)+"\t"+str((pos+1)*2000))


# Contact peaks file by the plasmid of the specie:

# file="/media/axel/EVO/dicty_lachancea_analysis_202408/Fermentati_peaks.bed"
# file="/media/axel/EVO/dicty_lachancea_analysis_202408/dicty_peaks.bed"
file="/media/axel/EVO/dicty_lachancea_analysis_202408/waltii_peaks.bed"

df= pd.read_table(file)

#rename of columns:
df.rename(columns = {'chrom1':'Chromosome'}, inplace = True)
df.rename(columns = {'start1':'Start'}, inplace = True)
df.rename(columns = {'end1':'End'}, inplace = True)
  
    
# TY regions or sub-telos or long gnes :
# file2="/home/axel/Bureau/YEAST/TY/positions_TY.txt"
# file2="/home/axel/Bureau/YEAST/sacCer3.sub_telo.txt2.30000.txt2"

# file2='/media/axel/d8e9c572-02dc-435e-a7f2-879cd0c77c86/dicty_lachancea_analysis_202408/GCA_900074765.1_LAFE0_genomic.gff.long.genes.7kb'
# file2='/media/axel/d8e9c572-02dc-435e-a7f2-879cd0c77c86/dicty_lachancea_analysis_202408/GCF_000004695.1_dicty_2.7_genomic.gff.long.genes.10kb'
file2='/media/axel/EVO/Lachencea/genome/waltii/embl_waltii/waltii_long_gene_7kb.txt'
       
df2= pd.read_table(file2)

# dictionary of sizes of chromosomes
# sizes = pd.read_table("/home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/sacCer3.chr_sizes.txt",
#                       header=None, delimiter="\t")

# sizes = pd.read_table("/media/axel/d8e9c572-02dc-435e-a7f2-879cd0c77c86/Lachencea/genome/fermentati/genome_fermentati.chrom.sizes.txt",
#                       header=None, delimiter="\t")

sizes = pd.read_table("/media/axel/d8e9c572-02dc-435e-a7f2-879cd0c77c86/Lachencea/hicstuff/FG0635/FG0635_waltii.chr.txt2",
                      header=None, delimiter="\t")

zipbObj = zip(sizes[0], sizes[1])
sizes_dist = dict(zipbObj)

zipbObj = zip(sizes[0], sizes[2])
names_dict = dict(zipbObj)


# translation des noms de chrom
df2['Chromosome'] = df2['Chromosome'].apply(lambda x: names_dict[x])

# conversion into pyranges objects
p1=pr.PyRanges(df)
p2=pr.PyRanges(df2)

len(p1)
len(p2)

# computation of intersections: 
pi=p2.intersect(p1, how="first")
len(pi)

# random group of peaks 
sizes = np.array(np.abs(df['Start']-df['End']))
mean_size = int(np.mean(sizes))

list_random=[]
for r in range(1000):
    p_r = pr.random(n=len(p1), length=mean_size, 
                    chromsizes=sizes_dist, int64=False)
    pri=p_r.intersect(p2, how="first")
    list_random.append(len(pri))


len(list_random)
np.mean(list_random)

# plot:
plt.hist(list_random,range(np.max(list_random)+2))
plt.xlabel("Number of overlaps")
plt.ylabel("Number of occurences")
plt.title("1000 Random groups")

  
array_random=np.array(list_random)
len(array_random[array_random>=2])

