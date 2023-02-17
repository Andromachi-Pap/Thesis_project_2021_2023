# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 01:02:31 2022

@author: User
"""


#SCRIPT TO PROCESS DESEQ RESULT TABLES AND GAIN UP-REGULATED AND DOWN-REGULATED GENOMIC ELEMENT IDs
 

import os
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
#pip install matplotlib-venn


path = os.chdir("C:/Users/User/Desktop/Msc_Bioinformatics_ 2020-2022/Master Thesis 2021-2022/Comparisons_rnaseq")
path

#open one result files from DESeq separately and read them into pandas dataframes
res_vibrio_invitro =  pd.read_csv("C:/Users/User/Desktop/Msc_Bioinformatics_ 2020-2022/Master Thesis 2021-2022/nf_rnaseq_invitro/results_vibrio_invitro_nfcore.tsv", sep='\t', header=0)
print(res_vibrio_invitro)
res_vibrio_invitro.rename( columns={'Unnamed: 0':'gene_id'}, inplace=True )
print(res_vibrio_invitro.columns)
res_pic_invitro = pd.read_csv("C:/Users/User/Desktop/Msc_Bioinformatics_ 2020-2022/Master Thesis 2021-2022/nf_rnaseq_invitro/results_PIC_invitro_nfcore.tsv", sep='\t', header=0)
res_pic_invitro.rename( columns={'Unnamed: 0':'gene_id'}, inplace=True )
print(res_pic_invitro)


res_vibrio_invivo =  pd.read_csv("C:/Users/User/Desktop/Msc_Bioinformatics_ 2020-2022/Master Thesis 2021-2022/nf_rnaseq_invivo/results_vibrio_invivo_nfcore.tsv", sep='\t', header=0)
res_vibrio_invivo.rename( columns={'Unnamed: 0':'gene_id'}, inplace=True )
#print(res_vibrio_invivo)

res_pic_invivo = pd.read_csv("C:/Users/User/Desktop/Msc_Bioinformatics_ 2020-2022/Master Thesis 2021-2022/nf_rnaseq_invivo/results_PIC_invivo_nfcore.tsv", sep='\t', header=0)
res_pic_invivo.rename( columns={'Unnamed: 0':'gene_id'}, inplace=True )
#print(res_pic_invivo)



#Let's identify the statistically significant and differentially expressed elements
#Up-regulated genes
#Take the genes (rows of the matrices) that have a padj < 0.1 and log2FoldChange > 0
up_vibrio_invitro = res_vibrio_invitro[(res_vibrio_invitro['padj'] < 0.1) & (res_vibrio_invitro['log2FoldChange'] > 0)]
up_vibrio_invitro.columns = ['gene_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue','padj']
up_vibrio_invitro.to_csv('up_vibrio_invitro_atac.csv' , index = False)
print(up_vibrio_invitro)

up_PIC_invitro =  res_pic_invitro[(res_pic_invitro['padj'] < 0.1) & (res_pic_invitro['log2FoldChange'] > 0)]
up_PIC_invitro.columns =  ['gene_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue','padj']
up_PIC_invitro.to_csv('up_PIC_invitro_atac.csv' , index = False)
print(up_PIC_invitro)


up_vibrio_invivo =   res_vibrio_invivo[(res_vibrio_invivo['padj'] < 0.1) & (res_vibrio_invivo['log2FoldChange'] > 0)]
up_vibrio_invivo.columns =  ['gene_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue','padj']
up_vibrio_invivo.to_csv('up_vibrio_invivo_atac.csv' , index = False)
print(up_vibrio_invivo)


up_PIC_invivo =  res_pic_invivo[(res_pic_invivo['padj'] < 0.1) & (res_pic_invivo['log2FoldChange'] > 0)]
up_PIC_invivo.columns =  ['gene_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue','padj']
up_PIC_invivo.to_csv('up_PIC_invivo.csv' , index = False)
print(up_PIC_invivo)

#Let's create lists of the UP-regulated genes on every treatment
gene_list_invitro_PIC_up = list(up_PIC_invitro["gene_id"])
gene_list_invitro_vibrio_up = list(up_vibrio_invitro["gene_id"])

gene_list_invivo_PIC_up = list(up_PIC_invivo["gene_id"])
gene_list_invivo_vibrio_up = list(up_vibrio_invivo["gene_id"])




#Let's save these lists into text files
with open(r'gene_list_atac_invitro_pic_up.txt', 'w') as fttp:
    fttp.write('\n'.join(gene_list_invitro_PIC_up))
    
with open(r'gene_list_atac_invitro_vibrio_up.txt', 'w') as pz:
    pz.write('\n'.join(gene_list_invitro_vibrio_up))

with open(r'gene_list_atac_invivo_pic_up.txt', 'w') as f:
    f.write('\n'.join(gene_list_invivo_PIC_up))
    
with open(r'gene_list_atac_invivo_vibrio_up.txt', 'w') as p:
    p.write('\n'.join(gene_list_invivo_vibrio_up))



#Down-regulated genes
#Take the genes (rows of the matrices) that have a padj < 0.1 and log2FoldChange < 0
down_vibrio_invitro = res_vibrio_invitro[(res_vibrio_invitro['padj'] < 0.1) & (res_vibrio_invitro['log2FoldChange'] < 0)]
down_vibrio_invitro.columns = ['gene_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue','padj']
down_vibrio_invitro.to_csv('down_vibrio_invitro_atac.csv' , index = False)
#print(down_vibrio_invitro)

down_PIC_invitro =  res_pic_invitro[(res_pic_invitro['padj'] < 0.1) & (res_pic_invitro['log2FoldChange'] < 0)]
down_PIC_invitro.columns =  ['gene_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue','padj']
down_PIC_invitro.to_csv('down_pic_invitro_atac.csv' , index = False)
#print(down_pic_invitro)

down_vibrio_invivo =   res_vibrio_invivo[(res_vibrio_invivo['padj'] < 0.1) & (res_vibrio_invivo['log2FoldChange'] < 0)]
down_vibrio_invivo.columns =  ['gene_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue','padj']
down_vibrio_invivo.to_csv('down_vibrio_invivo_atac.csv' , index = False)
#print(down_vibrio_invivo)

down_PIC_invivo =  res_pic_invivo[(res_pic_invivo['padj'] < 0.1) & (res_pic_invivo['log2FoldChange'] < 0)]
down_PIC_invivo.columns =  ['gene_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue','padj']
down_PIC_invivo.to_csv('down_pic_invivo_atac.csv' , index = False)
#print(down_pic_invivo)


#Let's create lists of the DOWN-regulated genes on every treatment
gene_list_invitro_PIC_down = list(down_PIC_invitro["gene_id"])
gene_list_invitro_vibrio_down = list(down_vibrio_invitro["gene_id"])

gene_list_invivo_PIC_down = list(down_PIC_invivo["gene_id"])
gene_list_invivo_vibrio_down = list(down_vibrio_invivo["gene_id"])


#Save these lists into txt files
with open(r'gene_list__atac_invitro_pic_down.txt', 'w') as yp:
    yp.write('\n'.join(gene_list_invitro_PIC_down))
    
with open(r'gene_list_atac_invitro_vibrio_down.txt', 'w') as rer:
    rer.write('\n'.join(gene_list_invitro_vibrio_down))

with open(r'gene_list_atac_invivo_pic_down.txt', 'w') as wz:
    wz.write('\n'.join(gene_list_invivo_PIC_down))
    
with open(r'gene_list_atac_invivo_vibrio_down.txt', 'w') as pu:
    pu.write('\n'.join(gene_list_invivo_vibrio_down))




#Gene ids with statistical significance  - padj < 0.1
vibrio_invitro_padj =  res_vibrio_invitro[(res_vibrio_invitro['padj'] < 0.1)]
vibrio_invivo_padj = res_vibrio_invivo[(res_vibrio_invivo['padj'] < 0.1)]
print(vibrio_invitro_padj)
vibrio_invitro_padj.to_csv('vibrio_invitro_rna_padj.csv' , index = False)
vibrio_invivo_padj.to_csv('vibrio_invivo_rna_padj.csv' , index = False)


PIC_invitro_padj = res_pic_invitro[(res_pic_invitro['padj'] < 0.1)] #empty
PIC_invivo_padj  = res_pic_invivo[(res_pic_invivo['padj'] < 0.1)]
print(PIC_invitro_padj)
print(PIC_invivo_padj)
PIC_invitro_padj.to_csv('PIC_invitro_rna_padj.csv' , index = False)
PIC_invivo_padj.to_csv('PIC_invivo_rna_padj.csv' , index = False)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#COMPARISONS - NETWORKS
#Let's do some comparisons between treatments for each setup (in-vitro, in-vivo)
#CORE IMMUNE RESPONSE NETWORK
#UP-regulated PIC vs Vibrio in-vitro
up_regulated_invitro = pd.merge(up_vibrio_invitro,up_PIC_invitro, how='inner', on= "gene_id")
gene_list_invitro_up = list(up_regulated_invitro["gene_id"])


#UP-regulated PIC vs Vibrio in-vivo
up_regulated_invivo = pd.merge(up_vibrio_invivo,up_PIC_invivo, how='inner', on= "gene_id")
gene_list_invivo_up = list(up_regulated_invivo["gene_id"])


#save the files into csv format
up_regulated_invitro.to_csv('up_regulated_invitro_rna.csv' , index = False)
up_regulated_invivo.to_csv('up_regulated_invivo_rna.csv' , index = False)

print(gene_list_invitro_up)
print(len(gene_list_invitro_up))

print(gene_list_invivo_up)
print(len(gene_list_invivo_up))

#store the created gene lists into text files
with open(r'rna_gene_list_invitro_up.txt', 'w') as fp:
    fp.write('\n'.join(gene_list_invitro_up))
    
with open(r'rna_gene_list_invivo_up.txt', 'w') as fpz:
    fpz.write('\n'.join(gene_list_invivo_up))
    
    
#DOWN-regulated PIC vs Vibrio in-vitro
down_regulated_invitro = pd.merge(down_vibrio_invitro,down_PIC_invitro, how='inner', on= "gene_id")
gene_list_invitro_down = list(down_regulated_invitro["gene_id"])


#DOWN-regulated PIC vs Vibrio in-vivo
down_regulated_invivo = pd.merge(down_vibrio_invivo,down_PIC_invivo, how='inner', on= "gene_id")
gene_list_invivo_down = list(down_regulated_invivo["gene_id"])

#save the files into csv format
down_regulated_invitro.to_csv('down_regulated_invitro_rna.csv' , index = False)
down_regulated_invivo.to_csv('down_regulated_invivo_rna.csv' , index = False)


#save the created gene lists into text files
with open(r'rna_gene_list_invitro_down.txt', 'w') as fpt:
    fpt.write('\n'.join(gene_list_invitro_down))
    
with open(r'rna_gene_list_invivo_down.txt', 'w') as fpo:
    fpo.write('\n'.join(gene_list_invivo_down))


#---------------------------------------------------------------------------------------------------

#UP-regulated PIC vs PIC (in-vitro,in-vivo)
up_regulated_pic_vs_pic = pd.merge(up_PIC_invitro,up_PIC_invivo, how='inner', on= "gene_id")
gene_list_pic_vs_pic_up = list(up_regulated_pic_vs_pic["gene_id"])

#UP-regulated Vibrio vs Vibrio (in-vitro,in-vivo)
up_regulated_vibrio_vs_vibrio = pd.merge(up_vibrio_invitro,up_vibrio_invivo, how='inner', on= "gene_id")
gene_list_vibrio_vs_vibrio_up = list(up_regulated_vibrio_vs_vibrio["gene_id"])

#save the files into csv format
up_regulated_pic_vs_pic.to_csv('up_regulated_pic_pic_rna.csv' , index = False)
up_regulated_vibrio_vs_vibrio.to_csv('up_regulated_vibrio_vibrio_rna.csv' , index = False)



#store the created gene lists into text files
with open(r'rna_gene_list_pic_vs_pic_up.txt', 'w') as p:
    p.write('\n'.join(gene_list_pic_vs_pic_up))
    
with open(r'rna_gene_list_vibrio_vs_vibrio_up.txt', 'w') as pz:
    pz.write('\n'.join(gene_list_vibrio_vs_vibrio_up))
    



#DOWN-regulated PIC vs PIC (in-vitro,in-vivo)
down_regulated_pic_vs_pic = pd.merge(down_PIC_invitro,down_PIC_invivo, how='inner', on= "gene_id")
gene_list_pic_vs_pic_down = list(down_regulated_pic_vs_pic["gene_id"])

#DOWN-regulated Vibrio vs Vibrio (in-vitro,in-vivo)
down_regulated_vibrio_vs_vibrio = pd.merge(down_vibrio_invitro,down_vibrio_invivo, how='inner', on= "gene_id")
gene_list_vibrio_vs_vibrio_down = list(down_regulated_vibrio_vs_vibrio["gene_id"])

#save the files into csv format
down_regulated_pic_vs_pic.to_csv('down_regulated_pic_pic_rna.csv' , index = False)
down_regulated_vibrio_vs_vibrio.to_csv('down_regulated_vibrio_vibrio_rna.csv' , index = False)



#store the created gene lists into text files
with open(r'rna_gene_list_pic_vs_pic_down.txt', 'w') as d:
    d.write('\n'.join(gene_list_pic_vs_pic_down))
    
with open(r'rna_gene_list_vibrio_vs_vibrio_down.txt', 'w') as dz:
    dz.write('\n'.join(gene_list_vibrio_vs_vibrio_down))
    
    
#print outcomes
print(gene_list_invitro_down)
print(len(gene_list_invitro_down))

print(gene_list_invivo_down)
print(len(gene_list_invivo_down))

#PIC
print(up_PIC_invitro)
print(up_PIC_invivo)

#VIBRIO
print(up_vibrio_invitro)
print(up_vibrio_invivo)

#------------------------------------------------------------------------------------------------
##Let's do some comparisons up and down regulated elements by between treatments in both setups
#VIRAL VS BACTERIAL NETWORK
up_pic_down_vibrio_invitro = pd.merge(up_PIC_invitro,down_vibrio_invitro, how='inner', on= "gene_id")
up_pic_down_vibrio_invitro_list = list(up_pic_down_vibrio_invitro["gene_id"])
up_pic_down_vibrio_invivo = pd.merge(up_PIC_invivo,down_vibrio_invivo, how='inner', on= "gene_id") #EMPTY

down_pic_up_vibrio_invitro =  pd.merge(down_PIC_invitro,up_vibrio_invitro, how='inner', on= "gene_id") #EMPTY
down_pic_up_vibrio_invivo =  pd.merge(down_PIC_invivo,up_vibrio_invivo, how='inner', on= "gene_id")
down_pic_up_vibrio_invivo_list = list(down_pic_up_vibrio_invivo["gene_id"])

#store the outcomes in csv files
up_pic_down_vibrio_invitro.to_csv('up_pic_down_vibrio_invitro_rna.csv' , index = False)
up_pic_down_vibrio_invivo.to_csv('up_pic_down_vibrio_invivo_rna.csv' , index = False)

down_pic_up_vibrio_invitro.to_csv('down_pic_up_vibrio_invitro_rna.csv' , index = False)
down_pic_up_vibrio_invivo.to_csv('down_pic_up_vibrio_invivo_rna.csv' , index = False)


#store the created gene lists into text files
with open(r'up_pic_down_vibrio_invitro_list.txt', 'w') as i:
    i.write('\n'.join(up_pic_down_vibrio_invitro_list))
    
with open(r'down_pic_up_vibrio_invivo_list.txt', 'w') as iz:
    iz.write('\n'.join(down_pic_up_vibrio_invivo_list))
    
#------------------------------------------------------------------------------------------------
#EXTRA NETWORK - NOT BIOLOGICALLY SIGNIFICANT?
#UP-regulated in both treatments in-Vitro AND/OR DOWN in both treatments in-vivo (AND OPPOSITE) 

#ΚΑΝΤΟ ΚΑΤΑΝΟΗΣΕ ΤΟ ΚΑΛΑ 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#PLOTTING - VENN-DIAGRAMS
#Let's create a Venn-Diagram for the common genes between in-vitro and in-vivo experiments for each treatment


#VENN DIAGRAMS
#UP

#Common up-regulated genes between in-vitro and in-vivo
#in vitro
venn2([set(list(up_vibrio_invitro["gene_id"])), set(list(up_PIC_invitro["gene_id"]))], set_labels = ('Vibrio', 'PIC'), set_colors= ("sandybrown", "slategray"), alpha = 0.7)
plt.title("Common up-regulated genes between Vibrio and PIC treatments in-vitro")
plt.savefig('Vibrio_vs_PIC_invitro_up_rna.pdf')
plt.show()


#in vivo
venn2([set(list(up_vibrio_invivo["gene_id"])), set(list(up_PIC_invivo["gene_id"]))], set_labels = ('Vibrio', 'PIC'), set_colors= ("darkorange", "lightslategray"), alpha = 0.7)
plt.title("Common up-regulated genes between Vibrio and PIC treatments in-vivo")
plt.savefig('Vibrio_vs_PIC_invivo_up_rna.pdf')
plt.show()


#PIC common genes between in-vitro, in-vivo
venn2([set(list(up_PIC_invitro["gene_id"])), set(list(up_PIC_invivo["gene_id"]))], set_labels = ('in-vitro', 'in-vivo'), set_colors= ("sandybrown", "slategray"), alpha = 0.7 )
plt.title("Common up-regulated genes between in-vitro and in-vivo conditions for PIC treatment")
plt.savefig('PIC_vs_PIC_up_rna.pdf')
plt.show()


#Vibrio common genes between in-vitro, in-vivo
venn2([set(list(up_vibrio_invitro["gene_id"])), set(list(up_vibrio_invivo["gene_id"]))], set_labels = ('in-vitro', 'in-vivo'),set_colors= ("sandybrown", "slategray"), alpha = 0.7)
plt.title("Common up-regulated genes between in-vitro and in-vivo conditions for Vibrio treatment")
plt.savefig('Vibrio_vs_Vibrio_up_rna.pdf')
plt.show()



#DOWN
#in- vitro
#Common down-regulated genes between in-vitro and in-vivo
venn2([set(list(down_vibrio_invitro["gene_id"])), set(list(down_PIC_invitro["gene_id"]))],set_labels = ('Vibrio', 'PIC'),  set_colors= ("chocolate", "steelblue"), alpha = 0.7)
plt.title("Common down-regulated genes between Vibrio and PIC treatments in-vitro")
plt.savefig('Vibrio_vs_PIC_invitro_down_rna.pdf')
plt.show()


#in vivo
venn2([set(list(down_vibrio_invivo["gene_id"])), set(list(down_PIC_invivo["gene_id"]))], set_labels = ('Vibrio', 'PIC'), set_colors= ("darkorange", "lightslategray"), alpha = 0.7)
plt.title("Common up-regulated genes between Vibrio and PIC treatments in-vivo")
plt.savefig('Vibrio_vs_PIC_invivo_down_rna.pdf')
plt.show()



#PIC common genes between in-vitro, in-vivo
venn2([set(list(down_PIC_invitro["gene_id"])), set(list(down_PIC_invivo["gene_id"]))], set_labels = ('in-vitro', 'in-vivo'),set_colors= ("chocolate", "steelblue"), alpha = 0.7);
plt.title("Common down-regulated genes between in-vitro and in-vivo conditions for PIC treatment")
plt.savefig('PIC_vs_PIC_down_rna.pdf')
plt.show()



#Vibrio common genes between in-vitro, in-vivo
venn2([set(list(down_vibrio_invitro["gene_id"])), set(list(down_vibrio_invivo["gene_id"]))], set_labels = ('in-vitro', 'in-vivo'), set_colors= ("chocolate", "steelblue"), alpha = 0.7)
plt.title("Common down-regulated genes between in-vitro and in-vivo conditions for Vibrio treatment")
plt.savefig('Vibrio_vs_Vibrio_down_rna.pdf')
plt.show()





#############################################################################################

#############################################################################################
#GENE LIST OF CHRIS 

#read Chris's gene list 
with open('genes_chris.txt') as file_in:
    chris_genes = []
    for line in file_in:
        line = line.rstrip("\n")
        chris_genes.append(line)
print(chris_genes)
#select the genes  of the above list from my results dataframes

print(res_vibrio_invitro)


#Create the subsets of the original results_Dataframes of DESeq with the gene IDs of interest

Vibrio_invitro= pd.DataFrame()
PIC_invitro = pd.DataFrame()

Vibrio_invivo = pd.DataFrame()
PIC_invivo = pd.DataFrame()

for gene in chris_genes:
    df_temp = res_vibrio_invitro.loc[res_vibrio_invitro['gene_id'] == gene]
    df_temp1 = res_pic_invitro.loc[res_pic_invitro['gene_id'] == gene]
    df_temp2 = res_vibrio_invivo.loc[res_vibrio_invivo['gene_id'] == gene]
    df_temp3 = res_pic_invivo.loc[res_pic_invivo['gene_id'] == gene]
    Vibrio_invitro = Vibrio_invitro.append(df_temp, ignore_index = True )
    PIC_invitro = PIC_invitro.append(df_temp1, ignore_index = True )
    Vibrio_invivo = Vibrio_invivo.append(df_temp2, ignore_index = True )
    PIC_invivo = PIC_invivo.append(df_temp3, ignore_index = True )
    
    
    
print(Vibrio_invitro)
print(PIC_invitro)
print(Vibrio_invivo)
print(PIC_invivo)



Vibrio_invitro.to_csv('Vibrio_invitro_chris.csv' , index = False)
PIC_invitro.to_csv('PIC_invitro_chris.csv' , index = False)    
Vibrio_invivo.to_csv('Vibrio_invivo_chris.csv' , index = False)
PIC_invivo.to_csv('PIC_invivo_chris.csv' , index = False)



#############################################################################################

#############################################################################################





