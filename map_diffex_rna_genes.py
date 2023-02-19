# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 00:11:26 2023

@author: User
"""



#ΜΑP THE RNA-SEQ GENE IDs WITH PADJ < 0.1 AND LOCATE THEM IN THE GTF FILE OF SPARUS AURATA
import os
import pandas as pd

path1 = os.chdir("Integrated_atacseq_rnaseq_folder")
path1

#We usesd the differentially expressed genes from the rna-seq analysis (padj < 0.1)
vibrio_invitro_rna_padj =  pd.read_csv("vibrio_invitro_rna_padj.csv", sep=',', header=0)
print(vibrio_invitro_rna_padj)

vibrio_invivo_rna_padj =  pd.read_csv("vibrio_invivo_rna_padj.csv", sep=',', header=0)
print(vibrio_invivo_rna_padj)

#empty df
pic_invitro_rna_padj =  pd.read_csv("pic_invitro_rna_padj.csv", sep=',', header=0)
print(pic_invitro_rna_padj)

pic_invivo_rna_padj =  pd.read_csv("pic_invivo_rna_padj.csv", sep=',', header=0)
print(pic_invivo_rna_padj)


#Let's extract the gene id's from the rna deseq result files
gene_list_invitro_PIC_rna = list(pic_invitro_rna_padj["gene_id"])
gene_list_invitro_vibrio_rna = list(vibrio_invitro_rna_padj["gene_id"])

gene_list_invivo_PIC_rna = list(pic_invivo_rna_padj["gene_id"])
gene_list_invivo_vibrio_rna = list(vibrio_invivo_rna_padj["gene_id"])




with open(r'gene_list_invitro_PIC_rna.txt', 'w', encoding='utf-8') as fpt:
    fpt.write('\n'.join(gene_list_invitro_PIC_rna))
    
with open(r'gene_list_invitro_vibrio_rna.txt', 'w', encoding='utf-8') as fpo:
    fpo.write('\n'.join(gene_list_invitro_vibrio_rna))

with open(r'gene_list_invivo_PIC_rna.txt', 'w', encoding='utf-8') as pt:
    pt.write('\n'.join(gene_list_invivo_PIC_rna))
    
with open(r'gene_list_invivo_vibrio_rna.txt', 'w', encoding='utf-8') as po:
    po.write('\n'.join(gene_list_invivo_vibrio_rna))




#####################################################################################################################

#We used a bed file of the genome, converted from the gtf
with open("Sparus_aurata.fSpaAur1.1.104.bed") as f:
        bed = list(f)
        
bed_sparus = [ line.split("\t") for line in bed ]


#Let's convert the whole bed file into a pandas df for better usage
bed_sparus_df = pd.DataFrame(bed_sparus)


#We only want specific columns of the initial bed to locate our genes
bed_to_parse = bed_sparus_df.iloc[:, 0:8]

#name the columns of the new df we created to locate our gene ids easily
bed_to_parse.columns= ["chr", "start", "end", "gene_id", "score", "strand", "source", "gene_type"]
print(bed_to_parse.iloc[0:13, :])


#we created a function that locates the gene ids into the bed dataframe made above for each of our gene lists
def map_rna_geneids(gene_list,bed_df):
    output_df = pd.DataFrame()
    for gene in gene_list:
        df_temp = bed_df.loc[bed_df["gene_id"] == gene]
        output_df = output_df.append(df_temp, ignore_index= True )
    return output_df
    

#let's locate all the gene ids from rna seq that were dif expr into our gtf file
PIC_invitro_rna_map = map_rna_geneids(gene_list_invitro_PIC_rna, bed_to_parse)
vibrio_invitro_rna_map = map_rna_geneids(gene_list_invitro_vibrio_rna, bed_to_parse)

PIC_invivo_rna_map = map_rna_geneids(gene_list_invivo_PIC_rna, bed_to_parse)
vibrio_invivo_rna_map = map_rna_geneids(gene_list_invivo_vibrio_rna, bed_to_parse)


#save the resulting dataframes into csv files
PIC_invitro_rna_map.to_csv('PIC_invitro_rna_map.csv' , index = False)
vibrio_invitro_rna_map.to_csv('vibrio_invitro_rna_map.csv' , index = False)

PIC_invivo_rna_map.to_csv('PIC_invivo_rna_map.csv' , index = False)
vibrio_invivo_rna_map.to_csv('vibrio_invivo_rna_map.csv' , index = False)
