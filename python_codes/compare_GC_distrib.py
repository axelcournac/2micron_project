#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 15:18:35 2021
@author: axel
To compare 
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')
import scipy.stats as st
import chromosight.kernels as ck

# data 
df1= pd.read_table('/home/axel/Bureau/GCcontent_HSC/gc_content_73HSC.txt', header=0)

df1_new = [pd.DataFrame(["HSC"] * len(df1["GC"])), 
           df1["GC"], 
           df1["GC-Skew"]]

headers = ["type", "GC", "GC-Skew"]
df1_new = pd.concat(df1_new, axis=1, keys=headers)
df1_new.columns=headers


df2= pd.read_table('/home/axel/Bureau/GCcontent_HSC/gc_stats.tsv', header=0)
df2["GC"]=df2["GC"]*100.0
df2["GCSKEW"]=df2["GCSKEW"]*100.0

df2_new = [pd.DataFrame(["whole_genome"] * len(df2["GC"])), 
           df2["GC"], 
           df2["GCSKEW"]]

headers = ["type", "GC", "GC-Skew"]
df2_new = pd.concat(df2_new, axis=1, keys=headers)
df2_new.columns=headers


frames = [df1_new, df2_new]
result = pd.concat(frames)

# Plot

sns.set_theme(style="white")
tips =  result

tips["all"] = ""
ax = sns.violinplot(x="all", y="GC-Skew", hue="type", data=tips, palette="muted", split=True)




