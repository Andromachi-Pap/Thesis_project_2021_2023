# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 23:23:40 2023

@author: User
"""

import os 
import pandas as pd
import numpy as np

###################### CREATE THE CIRCOS PLOT INPUT FILES ############################

#ADD LOG2FOLDCHANGE TO THE BED FILES OF EACH TREATMENT for RNA-seq
## SEPARATE THE BED FILES INTO UP AND DOWN REGULATED

path = os.chdir("Circos_plots_folder")
path

#read the atac and sorted rna bed files
rna_vibrio_invitro_bed = pd.read_csv("rna_vibrio_invitro_sorted.bed", sep = "\t", header = 0)
rna_vibrio_invivo_bed = pd.read_csv("rna_vibrio_invivo_sorted.bed", sep = "\t", header = 0)
rna_pic_invitro_bed = pd.read_csv("rna_pic_invitro_sorted.bed", sep = "\t", header = 0)
rna_pic_invivo_bed = pd.read_csv("rna_pic_invivo_sorted.bed", sep = "\t", header = 0)


#read the padj csv files of rna
rna_vibrio_invitro_padj = pd.read_csv("vibrio_invitro_rna_padj.csv", sep = ",", header = 0)
rna_vibrio_invivo_padj = pd.read_csv("vibrio_invivo_rna_padj.csv", sep = ",", header = 0)
rna_pic_invitro_padj = pd.read_csv("PIC_invitro_rna_padj.csv", sep = ",", header = 0)
rna_pic_invivo_padj = pd.read_csv("PIC_invivo_rna_padj.csv", sep = ",", header = 0)


#let's create smaller dataframes that only contain the GENE ID and the lfc
rna_vibrio_invitro_col = pd.concat([rna_vibrio_invitro_padj["gene_id"], rna_vibrio_invitro_padj["log2FoldChange"]], axis = 1)
rna_vibrio_invivo_col = pd.concat([rna_vibrio_invivo_padj["gene_id"], rna_vibrio_invivo_padj["log2FoldChange"]], axis = 1)
rna_pic_invitro_col = pd.concat([rna_pic_invitro_padj["gene_id"], rna_pic_invitro_padj["log2FoldChange"]], axis = 1)
rna_pic_invivo_col = pd.concat([rna_pic_invivo_padj["gene_id"], rna_pic_invivo_padj["log2FoldChange"]], axis = 1)


#add logfold change to bed files
rna_vibrio_invitro_lfc = pd.merge(rna_vibrio_invitro_bed,rna_vibrio_invitro_col, how = "inner", on="gene_id" )
rna_vibrio_invivo_lfc = pd.merge(rna_vibrio_invivo_bed,rna_vibrio_invivo_col, how = "inner", on="gene_id" )
rna_pic_invitro_lfc = pd.merge(rna_pic_invitro_bed,rna_pic_invitro_col, how = "inner", on="gene_id" )
rna_pic_invivo_lfc = pd.merge(rna_pic_invivo_bed,rna_pic_invivo_col, how = "inner", on="gene_id" )


#to make the files appropriate for circos input, we only keep specific columns chr, start,end, log2foldchange
rna_vibrio_invitro_comb = rna_vibrio_invitro_lfc[["chr", "start", "end", "log2FoldChange"]]
rna_vibrio_invivo_comb = rna_vibrio_invivo_lfc[["chr", "start", "end", "log2FoldChange"]]
rna_pic_invitro_comb = rna_pic_invitro_lfc[["chr", "start", "end", "log2FoldChange"]]
rna_pic_invivo_comb = rna_pic_invivo_lfc[["chr", "start", "end", "log2FoldChange"]]


#create a function that removes the rows that refer to scaffolds (not chromosomes) and add a chr preffix to each chromosome number
def prepare_circos_input(df):
    patternDel = "^CAAHF"
    filtr = df['chr'].str.contains(patternDel)
    new_df = df[~filtr]
    new_df['chr'] = 'chr' + new_df['chr'].astype(str)
    return new_df



rna_vibrio_invitro_bed_lfc = prepare_circos_input(rna_vibrio_invitro_comb)
rna_vibrio_invivo_bed_lfc = prepare_circos_input(rna_vibrio_invivo_comb)
rna_pic_invitro_bed_lfc = prepare_circos_input(rna_pic_invitro_comb)
rna_pic_invivo_bed_lfc = prepare_circos_input(rna_pic_invivo_comb)



#separate files into up and down-regulated

rna_vibrio_invitro_bedlfc_up = rna_vibrio_invitro_bed_lfc[(rna_vibrio_invitro_bed_lfc['log2FoldChange'] > 0)]
rna_vibrio_invivo_bedlfc_up = rna_vibrio_invivo_bed_lfc[(rna_vibrio_invivo_bed_lfc['log2FoldChange'] > 0)]
rna_pic_invitro_bedlfc_up = rna_pic_invitro_bed_lfc[(rna_pic_invitro_bed_lfc['log2FoldChange'] > 0)]
rna_pic_invivo_bedlfc_up = rna_pic_invivo_bed_lfc[(rna_pic_invivo_bed_lfc['log2FoldChange'] > 0)]


rna_vibrio_invitro_bedlfc_down = rna_vibrio_invitro_bed_lfc[(rna_vibrio_invitro_bed_lfc['log2FoldChange'] < 0)]
rna_vibrio_invivo_bedlfc_down = rna_vibrio_invivo_bed_lfc[(rna_vibrio_invivo_bed_lfc['log2FoldChange'] < 0)]
rna_pic_invitro_bedlfc_down = rna_pic_invitro_bed_lfc[(rna_pic_invitro_bed_lfc['log2FoldChange'] < 0)]
rna_pic_invivo_bedlfc_down = rna_pic_invivo_bed_lfc[(rna_pic_invivo_bed_lfc['log2FoldChange'] < 0)]



#save the files

rna_vibrio_invitro_bedlfc_up.to_csv("rna_vibrio_invitro_bedlfc_up.tsv",sep = "\t", header = 0, index = False)
rna_vibrio_invivo_bedlfc_up.to_csv("rna_vibrio_invivo_bedlfc_up.tsv",sep = "\t", header = 0, index = False)
rna_pic_invitro_bedlfc_up.to_csv("rna_pic_invitro_bedlfc_up.tsv", sep = "\t",header = 0, index = False)
rna_pic_invivo_bedlfc_up.to_csv("rna_pic_invivo_bedlfc_up.tsv", sep = "\t",header = 0, index = False)


rna_vibrio_invitro_bedlfc_down.to_csv("rna_vibrio_invitro_bedlfc_down.tsv", sep = "\t",header = 0, index = False)
rna_vibrio_invivo_bedlfc_down.to_csv("rna_vibrio_invivo_bedlfc_down.tsv", sep = "\t",header = 0, index = False)
rna_pic_invitro_bedlfc_down.to_csv("rna_pic_invitro_bedlfc_down.tsv", sep = "\t",header = 0, index = False)
rna_pic_invivo_bedlfc_down.to_csv("rna_pic_invivo_bedlfc_down.tsv",sep = "\t", header = 0, index = False)


#save the merged up-down dataframes into tsv files

rna_vibrio_invitro_bed_lfc.to_csv("rna_vibrio_invitro_padj_bedlfc.tsv", sep = "\t", header = 0, index = False)
rna_vibrio_invivo_bed_lfc.to_csv("rna_vibrio_invivo_padj_bedlfc.tsv", sep = "\t", header = 0, index = False)
rna_pic_invitro_bed_lfc.to_csv("rna_pic_invitro_padj_bedlfc.tsv", sep = "\t", header = 0, index = False)
rna_pic_invivo_bed_lfc.to_csv("rna_pic_invivo_padj_bedlfc.tsv", sep = "\t", header = 0, index = False)



#take the genome sizes file and make a karyotype input file for circos
genome_size = pd.read_csv("chrom_sizes_spaurata.genome", sep = "\t", header = None)
pattern_scaf = "^CAAHF"
filtr = genome_size.iloc[:,0].str.contains(pattern_scaf)
genome_size_df = genome_size[~filtr]
genome_size_df.iloc[:,0] = 'chr' + genome_size_df.iloc[:,0].astype(str)

karyotype_df = pd.DataFrame(np.random.randint(1,10, size=(24,7)))
karyotype_df.iloc[:,0] = "chr"
karyotype_df.iloc[:,1] = "-"
karyotype_df.iloc[:,2] = genome_size_df.iloc[:,0]
karyotype_df.iloc[:,3] = genome_size_df.iloc[:,0]
karyotype_df.iloc[:,4] = 0
karyotype_df.iloc[:,5] = genome_size_df.iloc[:,1]
karyotype_df.iloc[:,6] = ["vvdpurple","vdpurple","dpurple","purple","lpurple","vvdred","vdred","dred", "red", "lred","vvdblue","vdblue","dblue","blue","lblue","vvdgreen","vdgreen", "dgreen","green", "lgreen", "vvdyellow", "vdyellow", "dyellow", "yellow"]
print(karyotype_df)

#save the karyotype file
karyotype_df.to_csv("karyotype_file_spaurata.txt", sep = "\t", header = 0, index = False)


#find min - max of all files to define the value range of the circos plot
rna_vibrio_invitro_bedlfc_up["log2FoldChange"].max(axis=0)
rna_vibrio_invitro_bedlfc_up["log2FoldChange"].min(axis=0)
rna_vibrio_invitro_bedlfc_up.sort_values(by = "log2FoldChange").to_csv("rna_vibrio_invitro_up_sorted.tsv", sep = "\t",index = False)


rna_vibrio_invivo_bedlfc_up["log2FoldChange"].max(axis=0)
rna_vibrio_invivo_bedlfc_up["log2FoldChange"].min(axis=0)
rna_vibrio_invivo_bedlfc_up.sort_values(by = "log2FoldChange").to_csv("rna_vibrio_invivo_up_sorted.tsv", sep = "\t",index = False)

rna_pic_invivo_bedlfc_up["log2FoldChange"].max(axis=0)
rna_pic_invivo_bedlfc_up["log2FoldChange"].min(axis=0)
rna_pic_invivo_bedlfc_up.sort_values(by = "log2FoldChange").to_csv("rna_pic_invivo_up_sorted.tsv", sep = "\t",index = False)

rna_pic_invitro_bedlfc_up["log2FoldChange"].max(axis=0)
rna_pic_invitro_bedlfc_up["log2FoldChange"].min(axis=0)
rna_pic_invitro_bedlfc_up.sort_values(by = "log2FoldChange").to_csv("rna_pic_invitro_up_sorted.tsv", sep = "\t",index = False)


rna_vibrio_invitro_bedlfc_down["log2FoldChange"].max(axis=0)
rna_vibrio_invitro_bedlfc_down["log2FoldChange"].min(axis=0)
rna_vibrio_invitro_bedlfc_down.sort_values(by = "log2FoldChange").to_csv("rna_vibrio_invitro_down_sorted.tsv", sep = "\t",index = False)

rna_vibrio_invivo_bedlfc_down["log2FoldChange"].max(axis=0)
rna_vibrio_invivo_bedlfc_down["log2FoldChange"].min(axis=0)
rna_vibrio_invivo_bedlfc_down.sort_values(by = "log2FoldChange").to_csv("rna_vibrio_invivo_down_sorted.tsv", sep = "\t",index = False)

rna_pic_invitro_bedlfc_down["log2FoldChange"].max(axis=0)
rna_pic_invitro_bedlfc_down["log2FoldChange"].min(axis=0)
rna_pic_invitro_bedlfc_down.sort_values(by = "log2FoldChange").to_csv("rna_pic_invitro_down_sorted.tsv", sep = "\t",index = False)

rna_pic_invivo_bedlfc_down["log2FoldChange"].max(axis=0)
rna_pic_invivo_bedlfc_down["log2FoldChange"].min(axis=0)
rna_pic_invivo_bedlfc_down.sort_values(by = "log2FoldChange").to_csv("rna_pic_invivo_down_sorted.tsv", sep = "\t",index = False)



