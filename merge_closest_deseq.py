# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 01:15:52 2023

@author: User
"""

#MERGING ATAC-seq with RNA-seq
#PROCESS THE DATA TABLES FROM CLOSEST OF ATAC-RNA SEQ TO ADD LOG2FOLDCHANGE and find diffex genes

import os
import pandas as pd

path = os.chdir("ATAC_RNA_closest_processing_folder")
path

#Read the ATAC peak files
atac_vibrio_invitro =  pd.read_csv("vibrio_invitro_atac_padj.csv", sep=',', header=0)
#empty dataframe for pic in-vitro ATAC
atac_vibrio_invivo = pd.read_csv("vibrio_invivo_atac_padj.csv", sep=',', header=0)
atac_pic_invivo = pd.read_csv("PIC_invivo_atac_padj.csv", sep=',', header=0)


#Read the RNA gene ID files
rna_vibrio_invitro =  pd.read_csv("vibrio_invitro_rna_padj.csv", sep=',', header=0)
rna_vibrio_invivo =  pd.read_csv("vibrio_invivo_rna_padj.csv", sep=',', header=0)
#empty df for ATAC on PIC in-vitro, so we don't use the corresponding RNA-seq file
rna_pic_invivo = pd.read_csv("PIC_invivo_rna_padj.csv", sep=',', header=0)




#Read the bedtools-closest outcome of rnaseq elements mapped on atac peaks
closest_vibrio_invitro =  pd.read_csv("vibrio_invitro_closest_less_1Mb.bed", sep='\t', header = None)
closest_vibrio_invivo = pd.read_csv("vibrio_invivo_closest_less_1Mb.bed", sep='\t', header = None)
closest_pic_invivo = pd.read_csv("pic_invivo_closest_less_1Mb.bed", sep='\t', header = None)


#Name the columns of these dataframes for better handling
columns_vibrio_invitro = ['chr_rna', 'start_rna', 'end_rna', 'RNA_geneid', 'score1', 'strand', 'source', 'genetype', 'score2', 'chr_atac', 'start_atac', 'end_atac', 'gene_id', 'strand', 'length', 'distance' ]
columns_rest =  ["chr_rna", "start_rna", "end_rna", "RNA_geneid", "score1", "strand", "source", "genetype", "score2", "chr_atac", "start_atac", "end_atac", "gene_id", "strand",  "distance"  ]



closest_vibrio_invitro.set_axis(columns_vibrio_invitro, axis = 1 , inplace = True)
closest_vibrio_invivo.set_axis(columns_rest, axis =1, inplace = True)
closest_pic_invivo.set_axis(columns_rest, axis = 1, inplace = True)

#ADD LOGD2FOLDCHANGE TO THE DATAFRAMES
#for every interval ID that appears in the closest dataframes, the corresponding log2foldchange will be added as the last column 

#First, let's create smaller dataframes of the deseq atac data that only have the gene id and logfoldchange columns since we don't need the others
atac_vibrio_invitro_log = pd.concat([atac_vibrio_invitro.iloc[:,0],atac_vibrio_invitro.iloc[:, 2]], axis = 1)
atac_vibrio_invivo_log = pd.concat([atac_vibrio_invivo.iloc[:,0],atac_vibrio_invivo.iloc[:, 2]], axis = 1)
atac_pic_invivo_log = pd.concat([atac_pic_invivo.iloc[:,0],atac_pic_invivo.iloc[:, 2]], axis = 1)

#Let's do the same for the RNA seq gene IDs
rna_vibrio_invitro_log = pd.concat([rna_vibrio_invitro.iloc[:,0],rna_vibrio_invitro.iloc[:, 2]], axis = 1)
rna_vibrio_invivo_log = pd.concat([rna_vibrio_invivo.iloc[:,0],rna_vibrio_invivo.iloc[:, 2]], axis = 1)
rna_pic_invivo_log = pd.concat([rna_pic_invivo.iloc[:,0],rna_pic_invivo.iloc[:, 2]], axis = 1)
#rename the columns for better handling
rna_vibrio_invitro_log_1 = rna_vibrio_invitro_log.rename(columns={'gene_id': 'RNA_geneid', 'log2FoldChange': 'lfc_rna'})
rna_vibrio_invivo_log_1 = rna_vibrio_invivo_log.rename(columns={'gene_id': 'RNA_geneid', 'log2FoldChange': 'lfc_rna'})
rna_pic_invivo_log_1  = rna_pic_invivo_log.rename(columns={'gene_id': 'RNA_geneid', 'log2FoldChange': 'lfc_rna'})

#Let's match the atac peak files to their log2foldchnge values
merged_log_atac_vibrio_invitro = pd.merge(atac_vibrio_invitro_log,closest_vibrio_invitro, how='inner', on= "gene_id")
merged_log_atac_vibrio_invivo = pd.merge(atac_vibrio_invivo_log,closest_vibrio_invivo, how='inner', on= "gene_id")
merged_log_atac_pic_invivo = pd.merge(atac_pic_invivo_log,closest_pic_invivo, how='inner', on= "gene_id")

print(merged_log_atac_vibrio_invitro)
print(merged_log_atac_vibrio_invivo)
print(merged_log_atac_pic_invivo)

#Let's add now the log2foldchange values to the RNA gene IDs
atac_rna_lfc_vibrio_invitro = pd.merge(rna_vibrio_invitro_log_1,merged_log_atac_vibrio_invitro, how='inner', on= "RNA_geneid")
atac_rna_lfc_vibrio_invivo = pd.merge(rna_vibrio_invivo_log_1,merged_log_atac_vibrio_invivo, how='inner', on= "RNA_geneid")
atac_rna_lfc_pic_invivo = pd.merge(rna_pic_invivo_log_1,merged_log_atac_pic_invivo, how='inner', on= "RNA_geneid")

print(atac_rna_lfc_vibrio_invitro.columns)



#save the results into files
atac_rna_lfc_vibrio_invitro.to_csv("closest_vibrio_invitro_lfc.tsv", sep = "\t", index = False)
atac_rna_lfc_vibrio_invivo.to_csv("closest_vibrio_invivo_lfc.tsv", sep = "\t", index = False)
atac_rna_lfc_pic_invivo.to_csv("closest_pic_invivo_lfc.tsv", sep = "\t", index = False)




#UP AND DOWN-REGULATED DIFF ACCESSIBLE PEAKS
#UP
#Let's find the UP-regulated peaks from the files we created above
up_vibrio_invitro = atac_rna_lfc_vibrio_invitro[(atac_rna_lfc_vibrio_invitro['log2FoldChange'] > 0)]
up_vibrio_invivo = atac_rna_lfc_vibrio_invivo[(atac_rna_lfc_vibrio_invivo['log2FoldChange'] > 0) ]
up_pic_invivo = atac_rna_lfc_pic_invivo[(atac_rna_lfc_pic_invivo['log2FoldChange'] > 0)]
len(up_vibrio_invitro["gene_id"].unique())
len(up_vibrio_invivo["gene_id"].unique())
len(up_pic_invivo["gene_id"].unique())

#DOWN
#Let's find the down-regulated peaks from the files we created above
down_vibrio_invitro = atac_rna_lfc_vibrio_invitro[(atac_rna_lfc_vibrio_invitro['log2FoldChange'] < 0)]
down_vibrio_invivo = atac_rna_lfc_vibrio_invivo[(atac_rna_lfc_vibrio_invivo['log2FoldChange'] < 0)]
down_pic_invivo = atac_rna_lfc_pic_invivo[(atac_rna_lfc_pic_invivo['log2FoldChange'] < 0)]

#save the results into tsv files 
up_vibrio_invitro.to_csv("atac_vibrio_invitro_UP.tsv", index = False, sep = "\t")
up_vibrio_invivo.to_csv("atac_vibrio_invivo_UP.tsv", index = False,  sep = "\t")
up_pic_invivo.to_csv("atac_pic_invivo_UP.tsv", index = False, sep = "\t")

down_vibrio_invitro.to_csv("atac_vibrio_invitro_DOWN.tsv", index = False, sep = "\t")
down_vibrio_invivo.to_csv("atac_vibrio_invivo_DOWN.tsv", index = False, sep = "\t")
down_pic_invivo.to_csv("atac_pic_invivo_DOWN.tsv", index = False, sep = "\t")



#UP-regulated ATAC peaks and UP-regulated gene elements of RNA-seq
up_up_vibrio_invitro = atac_rna_lfc_vibrio_invitro[(atac_rna_lfc_vibrio_invitro['log2FoldChange'] > 0) & (atac_rna_lfc_vibrio_invitro['lfc_rna'] > 0)]
up_up_vibrio_invivo = atac_rna_lfc_vibrio_invivo[(atac_rna_lfc_vibrio_invivo['log2FoldChange'] > 0)  & (atac_rna_lfc_vibrio_invivo['lfc_rna'] > 0)]
up_up_pic_invivo = atac_rna_lfc_pic_invivo[(atac_rna_lfc_pic_invivo['log2FoldChange'] > 0) & (atac_rna_lfc_pic_invivo['lfc_rna'] > 0)]
print(up_up_vibrio_invitro)
print(up_up_vibrio_invivo)
print(up_up_pic_invivo)

#UP-regulated ATAC peaks and DOWN-regulated gene elements of RNA-seq
up_down_vibrio_invitro = atac_rna_lfc_vibrio_invitro[(atac_rna_lfc_vibrio_invitro['log2FoldChange'] > 0) & (atac_rna_lfc_vibrio_invitro['lfc_rna'] < 0)] 
up_down_vibrio_invivo = atac_rna_lfc_vibrio_invivo[(atac_rna_lfc_vibrio_invivo['log2FoldChange'] > 0) & (atac_rna_lfc_vibrio_invivo['lfc_rna'] < 0)]
up_down_pic_invivo = atac_rna_lfc_pic_invivo[(atac_rna_lfc_pic_invivo['log2FoldChange'] > 0) & (atac_rna_lfc_pic_invivo['lfc_rna'] < 0)]


#Save these records into files
up_up_vibrio_invitro.to_csv("atac_rna_vibrio_invitro_UP_UP.tsv", index = False, sep = "\t")
up_up_vibrio_invivo.to_csv("atac_rna_vibrio_invivo_UP_UP.tsv", index = False, sep = "\t")
up_up_pic_invivo.to_csv("atac_rna_pic_invivo_UP_UP.tsv", index = False, sep = "\t")

up_down_vibrio_invitro.to_csv("atac_rna_vibrio_invitro_UP_DOWN.tsv", index = False, sep = "\t")
up_down_vibrio_invivo.to_csv("atac_rna_vibrio_invivo_UP_DOWN.tsv", index = False, sep = "\t")
up_down_pic_invivo.to_csv("atac_rna_pic_invivo_UP_DOWN.tsv", index = False, sep = "\t")


#RNA GENE ID lists for UP-regulated ATAC peaks vs UP/DOWN-regulated gene elements of RNA-seq
up_up_vibrio_invitro_genes = up_up_vibrio_invitro["RNA_geneid"].unique()
up_up_vibrio_invivo_genes = up_up_vibrio_invivo["RNA_geneid"].unique()
up_up_pic_invivo_genes = up_up_pic_invivo["RNA_geneid"].unique()

up_down_vibrio_invitro_genes = up_down_vibrio_invitro["RNA_geneid"].unique()
up_down_vibrio_invivo_genes = up_down_vibrio_invivo["RNA_geneid"].unique()
up_down_pic_invivo_genes = up_down_pic_invivo["RNA_geneid"].unique()


#save these lists into text  files
with open(r'up_up_vibrio_invitro_genes.txt', 'w') as i:
    i.write('\n'.join(up_up_vibrio_invitro_genes))
 
with open(r'up_up_vibrio_invivo_genes.txt', 'w') as ia:
    ia.write('\n'.join(up_up_vibrio_invivo_genes))

with open(r'up_up_pic_invivo_genes.txt', 'w') as ib:
    ib.write('\n'.join(up_up_pic_invivo_genes))
    
with open(r'up_down_vibrio_invitro_genes.txt', 'w') as ic:
    ic.write('\n'.join(up_down_vibrio_invitro_genes))

with open(r'up_down_vibrio_invivo_genes.txt', 'w') as idd:
    idd.write('\n'.join(up_down_vibrio_invivo_genes))

with open(r'up_down_pic_invivo_genes.txt', 'w') as ie:
    ie.write('\n'.join(up_down_pic_invivo_genes))


    
#create files that have only the coordinates of the UP-regulated ATAC peaks and UP/DOWN-regulated gene elements
coordinates_up_up_vibrio_invitro = up_up_vibrio_invitro[['RNA_geneid', 'gene_id', 'chr_rna', 'start_rna', 'end_rna', 'chr_atac', 'start_atac', 'end_atac']]
coordinates_up_up_vibrio_invivo = up_up_vibrio_invivo[['RNA_geneid', 'gene_id', 'chr_rna', 'start_rna', 'end_rna', 'chr_atac', 'start_atac', 'end_atac']]
coordinates_up_up_pic_invivo = up_up_pic_invivo[['RNA_geneid', 'gene_id', 'chr_rna', 'start_rna', 'end_rna', 'chr_atac', 'start_atac', 'end_atac']]

coordinates_up_down_vibrio_invitro = up_down_vibrio_invitro[['RNA_geneid', 'gene_id', 'chr_rna', 'start_rna', 'end_rna', 'chr_atac', 'start_atac', 'end_atac']]
coordinates_up_down_vibrio_invivo = up_down_vibrio_invivo[['RNA_geneid', 'gene_id', 'chr_rna', 'start_rna', 'end_rna', 'chr_atac', 'start_atac', 'end_atac']]
coordinates_up_down_pic_invivo = up_down_pic_invivo[['RNA_geneid', 'gene_id', 'chr_rna', 'start_rna', 'end_rna', 'chr_atac', 'start_atac', 'end_atac']]

#save these coordinate into files
coordinates_up_up_vibrio_invitro.to_csv("atac_rna_vibrio_invitro_UP_UP_coord.tsv", index = False, sep = "\t")
coordinates_up_up_vibrio_invivo.to_csv("atac_rna_vibrio_invivo_UP_UP_coord.tsv", index = False, sep = "\t")
coordinates_up_up_pic_invivo.to_csv("atac_rna_pic_invivo_UP_UP_coord.tsv", index = False, sep = "\t")

coordinates_up_down_vibrio_invitro.to_csv("atac_rna_vibrio_invitro_UP_DOWN_coord.tsv", index = False, sep = "\t")
coordinates_up_down_vibrio_invivo.to_csv("atac_rna_vibrio_invivo_UP_DOWN_coord.tsv", index = False, sep = "\t")
coordinates_up_down_pic_invivo.to_csv("atac_rna_pic_invivo_UP_DOWN_coord.tsv", index = False, sep = "\t")




