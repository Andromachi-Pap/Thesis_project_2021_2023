# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 01:13:39 2023

@author: User
"""


import os
import pandas as pd


#SORT BED FILES FROM ATAC AND RNA SEQ FOR BEDTOOLS PROCESSES


#read BED Files
cwd = os.getcwd()
os.chdir("Working_directory")

#We show indicative file names that we used for this analysis
rna_vibrio_invitro = pd.read_csv("vibrio_invitro_rna_map.bed", header = 0, sep = ",")
rna_vibrio_invivo = pd.read_csv("vibrio_invivo_rna_map.bed", header = 0, sep = ",")

rna_pic_invitro = pd.read_csv("PIC_invitro_rna_map.bed", header = 0, sep = ",")
rna_pic_invivo = pd.read_csv("PIC_invivo_rna_map.bed", header = 0, sep = ",")


atac_vibrio_invitro = pd.read_csv("Vibrio_invitro.bed", header = 0, sep = ",")
#slice the dataframe to take the columns needed
atac_vibrio_invitro_bed = atac_vibrio_invitro.iloc[:,0:6]
#rename the Chromosome column from Chr to chr for easier parsing afterwards
atac_vibrio_invitro_bed.rename(columns={'Chr':'chr'}, inplace=True)

atac_vibrio_invivo = pd.read_csv("Vibrio_invivo.bed", header = 0, sep = ",")
atac_vibrio_invivo.rename(columns={'Chr':'chr'}, inplace=True)

#atac_pic_invitro --- EMPTY 

atac_pic_invivo = pd.read_csv("PIC_invivo.bed", header = 0, sep = ",")
atac_pic_invivo.rename(columns={'Chr':'chr'}, inplace=True)



#Problem with sorting the pandas dataframes because the list of chromosomes is parsed as string, containing float and text
#Let's fix that by converting the strings to numbers and keep the text-chromosomes apart to sort correctly and merge them again

def sort_chr(df):
    #find the rows that contain chromosomes with text
    char = df[df.chr.str.contains(r'[CA]')]
    #keep their indexes
    list_char = char.index
    #remove these rows from the whole dataframe
    num = df.drop(axis =0, index = list_char)
    #convert the chromosome column to numbers
    num["chr"] = pd.to_numeric(num["chr"])
    #sort the chromosomes by number
    sorted_chr = num.sort_values("chr")
    #merge the character-chromosomes with the numeric ones to get the whole dataframe again
    total_df = pd.concat([sorted_chr, char], axis = 0)
    return total_df


##Sort bed files for the betools-closest process
rna_vibrio_invitro_sorted = sort_chr(rna_vibrio_invitro)
rna_vibrio_invivo_sorted = sort_chr(rna_vibrio_invivo)

rna_pic_invitro_sorted = sort_chr(rna_pic_invitro)
rna_pic_invivo_sorted = sort_chr(rna_pic_invivo)

atac_vibrio_invitro_sorted = sort_chr(atac_vibrio_invitro_bed)


atac_vibrio_invivo_sorted = sort_chr(atac_vibrio_invivo)
atac_pic_invivo_sorted = sort_chr(atac_pic_invivo)

list(atac_vibrio_invivo_sorted["chr"])

list(atac_vibrio_invivo_sorted["chr"])
#to make the atac peaks files readable by bedtools, we need the chr on the first column and their coordinates on 2-3 columns


def move_column_inplace(df):
    gene_id = df["Geneid"]
    col = df.pop("Geneid")
    df.insert(3, "Geneid", gene_id)
    return df

atac_vibrio_invitro_sorted_ = move_column_inplace(atac_vibrio_invitro_sorted)
atac_vibrio_invivo_sorted_ = move_column_inplace(atac_vibrio_invivo_sorted)

atac_pic_invivo_sorted_ = move_column_inplace(atac_pic_invivo_sorted)


#let's make all the column names identical to the ones of the rna experiments for easier parsing 
atac_vibrio_invitro_sorted_.columns = atac_vibrio_invitro_sorted_.columns.str.lower()
atac_vibrio_invivo_sorted_.columns = atac_vibrio_invivo_sorted_.columns.str.lower()

atac_pic_invivo_sorted_.columns = atac_pic_invivo_sorted_.columns.str.lower()


#We need to make sure that the columns of chromosome coordinates are sorted within each chromosome

def sort_within_chromosome(df):
    output_df = pd.DataFrame()
    #group smaller slices of the initial dataframe for each chromosome 
    grouped = df.groupby(df.chr)
    #for each sub-group of the initial df that is about one chromosome, sort the coordinates and append the chromosome group to a temp df and then to the output df
    for g in grouped.groups:
        df_temp = grouped.get_group(g).sort_values("start")
        output_df = output_df.append(df_temp, ignore_index= True )
    return output_df

rna_vibrio_invitro_sorted_final = sort_within_chromosome(rna_vibrio_invitro_sorted)
rna_vibrio_invivo_sorted_final = sort_within_chromosome(rna_vibrio_invivo_sorted)

rna_pic_invitro_sorted_final = sort_within_chromosome(rna_pic_invitro_sorted)
rna_pic_invivo_sorted_final = sort_within_chromosome(rna_pic_invivo_sorted)



atac_vibrio_invitro_sorted_final = sort_within_chromosome(atac_vibrio_invitro_sorted_)
atac_vibrio_invivo_sorted_final = sort_within_chromosome(atac_vibrio_invivo_sorted_)

atac_pic_invivo_sorted_final = sort_within_chromosome(atac_pic_invivo_sorted_)



#Save the files
rna_vibrio_invitro_sorted_final.to_csv("rna_vibrio_invitro_sorted.bed", sep = "\t", index = False)
rna_vibrio_invivo_sorted_final.to_csv("rna_vibrio_invivo_sorted.bed", sep = "\t", index = False)
rna_pic_invitro_sorted_final.to_csv("rna_pic_invitro_sorted.bed", sep = "\t", index = False)
rna_pic_invivo_sorted_final.to_csv("rna_pic_invivo_sorted.bed", sep = "\t", index = False)

atac_vibrio_invitro_sorted_final.to_csv("atac_vibrio_invitro_sorted.bed", sep = "\t", index = False)
atac_vibrio_invivo_sorted_final.to_csv("atac_vibrio_invivo_sorted.bed", sep ="\t", index = False)
atac_pic_invivo_sorted_final.to_csv("atac_pic_invivo_sorted.bed", sep ="\t", index = False)



#Let's match the choromosomes that appear on the atac peaks with the rna genes
#for every chromosome that appears in atac peaks file, the same number of chromosome will be taken from the next treatment of rna-seq
chrom_vibrio_invitro_atac = list(atac_vibrio_invitro_sorted_final.chr.unique())
chrom_vibrio_invivo_atac = list(atac_vibrio_invivo_sorted_final.chr.unique())
chrom_pic_invivo_atac = list(atac_pic_invivo_sorted_final.chr.unique())

#We created a function that does the following:
#for every chromosome that appears in atac peaks file, the same number of chromosome will be taken from the next treatment of rna-seq
def match_rna_atac(chr_list,df):
    output_df = pd.DataFrame()
    for chrom in chr_list:
        df_temp = df.loc[df["chr"] == chrom]
        output_df = output_df.append(df_temp, ignore_index= True )
    return output_df


rna_vibrio_invitro_matched = match_rna_atac(chrom_vibrio_invitro_atac, rna_vibrio_invitro_sorted_final)
rna_vibrio_invivo_matched = match_rna_atac(chrom_vibrio_invivo_atac, rna_vibrio_invivo_sorted_final)
rna_pic_invivo_matched = match_rna_atac(chrom_pic_invivo_atac, rna_pic_invivo_sorted_final)



#save the results into bed files
rna_vibrio_invitro_matched.to_csv("rna_vibrio_invitro_matched.bed", sep ="\t", index = False)
rna_vibrio_invivo_matched.to_csv("rna_vibrio_invivo_matched.bed", sep ="\t", index = False)
rna_pic_invivo_matched.to_csv("rna_pic_invivo_matched.bed", sep ="\t", index = False)


#REVERSE: SEARCH THE CHROMOSOME NUMBERS OF RNA-SEQ ON ATAC-SEQ
#Let's do the reverse process now: match the chromosome numbers of rna-seq on atac-seq files
chrom_vibrio_invitro_rna = list(rna_vibrio_invitro_sorted_final.chr.unique())
chrom_vibrio_invivo_rna = list(rna_vibrio_invivo_sorted_final.chr.unique())
chrom_pic_invivo_rna = list(rna_pic_invivo_sorted_final.chr.unique())

atac_vibrio_invitro_matched = match_rna_atac(chrom_vibrio_invitro_rna, atac_vibrio_invitro_sorted_final)
atac_vibrio_invivo_matched = match_rna_atac(chrom_vibrio_invivo_rna, atac_vibrio_invivo_sorted_final)
atac_pic_invivo_matched = match_rna_atac(chrom_pic_invivo_rna, atac_pic_invivo_sorted_final)



#save the results into bed files
atac_vibrio_invitro_matched.to_csv("atac_vibrio_invitro_matched.bed", sep ="\t", index = False)
atac_vibrio_invivo_matched.to_csv("atac_vibrio_invivo_matched.bed", sep ="\t", index = False)
atac_pic_invivo_matched.to_csv("atac_pic_invivo_matched.bed", sep ="\t", index = False)