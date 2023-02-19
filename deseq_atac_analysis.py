# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 01:02:31 2022

@author: User
"""


#SCRIPT TO PROCESS DESEQ RESULT TABLES AND GAIN UP-REGULATED AND DOWN-REGULATED PEAK FEATURE IDs FOR ATAC-SEQ DATA
 

import os
import pandas as pd


#Define the working directory
path = os.chdir("C:/Users/User/Desktop/atac_seq_deseq_folder")
path


#We gave indicative file names that represent each treatment dataset - user may name them differently
#open one-by-one the result files from DESeq separately and read them into pandas dataframes

res_vibrio_invitro =  pd.read_csv("results_vibrio_invitro_atac.tsv", sep='\t', header=0)
print(res_vibrio_invitro)
#rename the peak_id column because in the original dataframe there is no column name
res_vibrio_invitro.rename( columns={'Unnamed: 0':'peak_id'}, inplace=True )
print(res_vibrio_invitro.columns)
res_pic_invitro = pd.read_csv("results_PIC_invitro_atac.tsv", sep='\t', header=0)
#rename the peak_id column because in the original dataframe there is no column name
res_pic_invitro.rename( columns={'Unnamed: 0':'peak_id'}, inplace=True )
print(res_pic_invitro)


res_vibrio_invivo =  pd.read_csv("results_vibrio_invivo_atac.tsv", sep='\t', header=0)
#rename the peak_id column because in the original dataframe there is no column name
res_vibrio_invivo.rename( columns={'Unnamed: 0':'peak_id'}, inplace=True )
print(res_vibrio_invivo)
res_pic_invivo = pd.read_csv("results_PIC_invivo_atac.tsv", sep='\t', header=0)
#rename the peak_id column because in the original dataframe there is no column name
res_pic_invivo.rename( columns={'Unnamed: 0':'peak_id'}, inplace=True )
print(res_pic_invivo)


#Let's identify the statistically significant and differentially accessible elements

#Peak ids with statistical significance  - padj < 0.1
vibrio_invitro_padj =  res_vibrio_invitro[(res_vibrio_invitro['padj'] < 0.1)]
vibrio_invivo_padj = res_vibrio_invivo[(res_vibrio_invivo['padj'] < 0.1)]
print(vibrio_invitro_padj)
vibrio_invitro_padj.to_csv('vibrio_invitro_atac_padj.csv' , index = False)
vibrio_invivo_padj.to_csv('vibrio_invivo_atac_padj.csv' , index = False)
print(vibrio_invivo_padj)

PIC_invitro_padj = res_pic_invitro[(res_pic_invitro['padj'] < 0.1)] #empty for this dataset 
PIC_invivo_padj  = res_pic_invivo[(res_pic_invivo['padj'] < 0.1)]
print(PIC_invitro_padj)
print(PIC_invivo_padj)
PIC_invitro_padj.to_csv('PIC_invitro_atac_padj.csv' , index = False)
PIC_invivo_padj.to_csv('PIC_invivo_atac_padj.csv' , index = False)

#------------------------------------------------------------------------------------------------------------------------

#Differential accessibility analysis
#Up-regulated peaks
#Take the peaks (rows of the matrices) that have a padj < 0.1 and log2FoldChange > 0
up_vibrio_invitro = res_vibrio_invitro[(res_vibrio_invitro['padj'] < 0.1) & (res_vibrio_invitro['log2FoldChange'] > 0)]
up_vibrio_invitro.columns = ['peak_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue','padj']
up_vibrio_invitro.to_csv('up_vibrio_invitro_atac.csv' , index = False)
print(up_vibrio_invitro)

up_PIC_invitro =  res_pic_invitro[(res_pic_invitro['padj'] < 0.1) & (res_pic_invitro['log2FoldChange'] > 0)]
up_PIC_invitro.columns =  ['peak_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue','padj']
up_PIC_invitro.to_csv('up_PIC_invitro_atac.csv' , index = False)
print(up_PIC_invitro)


up_vibrio_invivo =   res_vibrio_invivo[(res_vibrio_invivo['padj'] < 0.1) & (res_vibrio_invivo['log2FoldChange'] > 0)]
up_vibrio_invivo.columns =  ['peak_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue','padj']
up_vibrio_invivo.to_csv('up_vibrio_invivo_atac.csv' , index = False)
print(up_vibrio_invivo)


up_PIC_invivo =  res_pic_invivo[(res_pic_invivo['padj'] < 0.1) & (res_pic_invivo['log2FoldChange'] > 0)]
up_PIC_invivo.columns =  ['peak_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue','padj']
up_PIC_invivo.to_csv('up_PIC_invivo_atac.csv' , index = False)
print(up_PIC_invivo)

#Let's create lists of the UP-regulated peaks on every treatment
peak_list_invitro_PIC_up = list(up_PIC_invitro["peak_id"])
peak_list_invitro_vibrio_up = list(up_vibrio_invitro["peak_id"])

peak_list_invivo_PIC_up = list(up_PIC_invivo["peak_id"])
peak_list_invivo_vibrio_up = list(up_vibrio_invivo["peak_id"])




#Let's save these lists into text files
with open(r'peak_list_atac_invitro_pic_up.txt', 'w') as fttp:
    fttp.write('\n'.join(peak_list_invitro_PIC_up))
    
with open(r'peak_list_atac_invitro_vibrio_up.txt', 'w') as pz:
    pz.write('\n'.join(peak_list_invitro_vibrio_up))

with open(r'peak_list_atac_invivo_pic_up.txt', 'w') as f:
    f.write('\n'.join(peak_list_invivo_PIC_up))
    
with open(r'peak_list_atac_invivo_vibrio_up.txt', 'w') as p:
    p.write('\n'.join(peak_list_invivo_vibrio_up))



#Down-regulated peaks
#Take the peaks (rows of the matrices) that have a padj < 0.1 and log2FoldChange < 0
down_vibrio_invitro = res_vibrio_invitro[(res_vibrio_invitro['padj'] < 0.1) & (res_vibrio_invitro['log2FoldChange'] < 0)]
down_vibrio_invitro.columns = ['peak_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue','padj']
down_vibrio_invitro.to_csv('down_vibrio_invitro_atac.csv' , index = False)
print(down_vibrio_invitro)

down_PIC_invitro =  res_pic_invitro[(res_pic_invitro['padj'] < 0.1) & (res_pic_invitro['log2FoldChange'] < 0)]
down_PIC_invitro.columns =  ['peak_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue','padj']
down_PIC_invitro.to_csv('down_pic_invitro_atac.csv' , index = False)
print(down_PIC_invitro)

down_vibrio_invivo =   res_vibrio_invivo[(res_vibrio_invivo['padj'] < 0.1) & (res_vibrio_invivo['log2FoldChange'] < 0)]
down_vibrio_invivo.columns =  ['peak_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue','padj']
down_vibrio_invivo.to_csv('down_vibrio_invivo_atac.csv' , index = False)
print(down_vibrio_invivo)

down_PIC_invivo =  res_pic_invivo[(res_pic_invivo['padj'] < 0.1) & (res_pic_invivo['log2FoldChange'] < 0)]
down_PIC_invivo.columns =  ['peak_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue','padj']
down_PIC_invivo.to_csv('down_pic_invivo_atac.csv' , index = False)
print(down_PIC_invivo)


#Let's create lists of the DOWN-regulated peaks for every treatment
peak_list_invitro_PIC_down = list(down_PIC_invitro["peak_id"])
peak_list_invitro_vibrio_down = list(down_vibrio_invitro["peak_id"])

peak_list_invivo_PIC_down = list(down_PIC_invivo["peak_id"])
peak_list_invivo_vibrio_down = list(down_vibrio_invivo["peak_id"])


#Save these lists into txt files
with open(r'peak_list_atac_invitro_pic_down.txt', 'w') as yp:
    yp.write('\n'.join(peak_list_invitro_PIC_down))
    
with open(r'peak_list_atac_invitro_vibrio_down.txt', 'w') as rer:
    rer.write('\n'.join(peak_list_invitro_vibrio_down))

with open(r'peak_list_atac_invivo_pic_down.txt', 'w') as wz:
    wz.write('\n'.join(peak_list_invivo_PIC_down))
    
with open(r'peak_list_atac_invivo_vibrio_down.txt', 'w') as pu:
    pu.write('\n'.join(peak_list_invivo_vibrio_down))




