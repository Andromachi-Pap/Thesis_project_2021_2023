# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 01:35:42 2023

@author: User
"""

#FIND ENHANCERS/SILENCERS FOR UP-REGULATED GENES
import os
import pandas as pd

path = os.chdir("enhancers_silencers_for_genes_folder")
path

#ENHANCERS
#Let's find out which up-regulated genes from RNA-seq are next to up-regulated atac peaks

#load the PIC-Vibrio up-regulated gene tables 
#in-vitro
with open("rna_gene_list_invitro_up.txt") as file_in_1:
    gene_list_invitro_up = []
    for line in file_in_1:
        line = line.rstrip("\n")
        gene_list_invitro_up.append(line)
print(gene_list_invitro_up)
print(len(gene_list_invitro_up))

#in-vivo
with open("rna_gene_list_invivo_up.txt") as file_in_2:
    gene_list_invivo_up = []
    for line in file_in_2:
        line = line.rstrip("\n")
        gene_list_invivo_up.append(line)
print(gene_list_invivo_up)
print(len(gene_list_invivo_up))


#load the PIC-Vibrio down-regulated gene tables 
#in-vitro
with open("rna_gene_list_invitro_down.txt") as file_in_3:
    gene_list_invitro_down = []
    for line in file_in_3:
        line = line.rstrip("\n")
        gene_list_invitro_down.append(line)
print(gene_list_invitro_down)
print(len(gene_list_invitro_down))

#in-vivo
with open("rna_gene_list_invivo_down.txt") as file_in_4:
    gene_list_invivo_down = []
    for line in file_in_4:
        line = line.rstrip("\n")
        gene_list_invivo_down.append(line)
print(gene_list_invivo_down)
print(len(gene_list_invivo_down))


#load the UP-regulated atac peak files
up_vibrio_invitro_atac = pd.read_table("atac_vibrio_invitro_UP.tsv")
up_vibrio_invivo_atac = pd.read_table("atac_vibrio_invivo_UP.tsv")
up_pic_invivo_atac = pd.read_table("atac_pic_invivo_UP.tsv")


print(len(up_vibrio_invitro_atac["gene_id"].unique()))
print(len(up_vibrio_invivo_atac["gene_id"].unique()))
print(len(up_pic_invivo_atac["gene_id"].unique()))


#Let's concatenate all the in-vivo tables together to test in-vitro and in-vivo elements separately
up_invivo_atac = pd.concat([up_vibrio_invivo_atac, up_pic_invivo_atac])

enhancers_invitro = pd.DataFrame()
enhancers_invivo = pd.DataFrame()

for gene in gene_list_invitro_up:
    df_temp = up_vibrio_invitro_atac.loc[up_vibrio_invitro_atac['RNA_geneid'] == gene]
    enhancers_invitro = enhancers_invitro.append(df_temp, ignore_index = True )

for gene_2 in gene_list_invivo_up:
    df_temp2 = up_invivo_atac.loc[up_invivo_atac["RNA_geneid"] == gene_2]
    enhancers_invivo = enhancers_invivo.append(df_temp2, ignore_index = True)
    
#let's save the results into tsv files
enhancers_invitro.to_csv("enhancers_invitro.tsv", sep = "\t", index = False)
enhancers_invivo.to_csv("enhancers_invivo.tsv", sep = "\t", index = False)


enhanced_genes_invitro = enhancers_invitro["RNA_geneid"].unique()
enhanced_genes_invivo = enhancers_invivo["RNA_geneid"].unique()

#store the created gene lists into text files
with open(r'rna_genes_list_enhancers_invitro.txt', 'w') as fp:
    fp.write('\n'.join(enhanced_genes_invitro))
    
with open(r'rna_genes_list_enhancers_invivo.txt', 'w') as fpz:
    fpz.write('\n'.join(enhanced_genes_invivo))



#SILENCERS type 1
#Let's find out which up-regulated genes from RNA-seq are next to down-regulated atac peaks

#load the DOWN-regulated atac peaks
#Vibrio in-vitro has 0 DOWN-regulated atac peaks
down_vibrio_invivo_atac = pd.read_table("atac_vibrio_invivo_DOWN.tsv")
down_pic_invivo_atac = pd.read_table("atac_pic_invivo_DOWN.tsv")

#let's concatenate the two treatments for easier handling of data in-vivo all together
down_invivo_atac = pd.concat([down_vibrio_invivo_atac, down_pic_invivo_atac])

silencers_invivo = pd.DataFrame()

for gene_3 in gene_list_invivo_up:
     df_temp = down_invivo_atac.loc[down_invivo_atac['RNA_geneid'] == gene_3]
     silencers_invivo = silencers_invivo.append(df_temp, ignore_index = True )
     
print(silencers_invivo)
print(silencers_invivo["RNA_geneid"].unique())

silenced_genes_invivo = silencers_invivo["RNA_geneid"].unique()

#store the created gene lists into text files
with open(r'rna_genes_list_silencers_invivo.txt', 'w') as fz:
    fz.write('\n'.join(silenced_genes_invivo))

#let's save the results into tsv files
silencers_invivo.to_csv("silencers_invivo.tsv", sep = "\t", index = False)



#SILENCERS type 2
##Let's find out which down-regulated genes from RNA-seq are next to up-regulated atac peaks
silencers2_invitro = pd.DataFrame()
silencers2_invivo = pd.DataFrame()

for gene in gene_list_invitro_down:
    df_temp = up_vibrio_invitro_atac.loc[up_vibrio_invitro_atac['RNA_geneid'] == gene]
    silencers2_invitro = silencers2_invitro.append(df_temp, ignore_index = True )

for gene_2 in gene_list_invivo_down:
    df_temp2 = up_invivo_atac.loc[up_invivo_atac["RNA_geneid"] == gene_2]
    silencers2_invivo = silencers2_invivo.append(df_temp2, ignore_index = True)
    
#let's save the results into tsv files
silencers2_invitro.to_csv("silencers2_invitro.tsv", sep = "\t", index = False)
silencers2_invivo.to_csv("silencers2_invivo.tsv", sep = "\t", index = False)


silenced2_genes_invitro = silencers2_invitro["RNA_geneid"].unique()
silenced2_genes_invivo = silencers2_invivo["RNA_geneid"].unique()

#store the created gene lists into text files
with open(r'rna_genes_list_silencers2_invitro.txt', 'w') as fp:
    fp.write('\n'.join(silenced2_genes_invitro))
    
with open(r'rna_genes_list_silencers2_invivo.txt', 'w') as fpz:
    fpz.write('\n'.join(silenced2_genes_invivo))

