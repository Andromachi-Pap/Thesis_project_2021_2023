# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 00:05:44 2023

@author: User
"""


#SCRIPT THAT MAPS THE DIFFERENTIALLY ACCESSIBLE PEAKS ONTO THE GENOME (ATAC-SEQ)

import os
import pandas as pd


#Load the differentially accessible peaks
#we used the ones we defined from deseq from nf-core run of atac-seq (padj < 0.1)
path = os.chdir("Working directory")
path

vibrio_invitro_padj =  pd.read_csv("vibrio_invitro_atac_padj.csv", sep=',', header=0)
print(vibrio_invitro_padj)

vibrio_invivo_padj =  pd.read_csv("vibrio_invivo_atac_padj.csv", sep=',', header=0)
print(vibrio_invivo_padj)

#empty df
pic_invitro_padj =  pd.read_csv("pic_invitro_atac_padj.csv", sep=',', header=0)
print(pic_invitro_padj)

pic_invivo_padj =  pd.read_csv("pic_invivo_atac_padj.csv", sep=',', header=0)
print(pic_invivo_padj)


#Map them into the genome using the consensus peaks result txt files
#Load the consensus peaks files - indicative names of the files we used
consensus_peaks_invitro =  pd.read_csv("consensus_peaks.mLb.clN.results_invitro.txt", sep='\t', header=0)
#print(consensus_peaks_invitro)

consensus_peaks_invivo = pd.read_csv("consensus_peaks.mLb.clN.results_invivo.txt", sep = "\t", header=0)
#print(consensus_peaks_invivo)

#For each of the interval ID's of every treatment, the coordinates of the intervals will be listed
vibrio_invitro_list = list(vibrio_invitro_padj["gene_id"])
vibrio_invivo_list = list(vibrio_invivo_padj["gene_id"])

pic_invivo_list = list(pic_invivo_padj["gene_id"])


#create new dataframes to insert each peak and its coordinates
Vibrio_invitro= pd.DataFrame()

Vibrio_invivo = pd.DataFrame()
PIC_invivo = pd.DataFrame()

for interval in vibrio_invitro_list:
     df_temp = consensus_peaks_invitro.loc[consensus_peaks_invitro['Geneid'] == interval]
     Vibrio_invitro = Vibrio_invitro.append(df_temp, ignore_index = True )

         
for interval in vibrio_invivo_list:
    df_temp1 = consensus_peaks_invivo.loc[consensus_peaks_invivo['Geneid'] == interval]
    Vibrio_invivo = Vibrio_invivo.append(df_temp1, ignore_index = True )
    
     
for interval in pic_invivo_list:
    df_temp2 = consensus_peaks_invivo.loc[consensus_peaks_invivo['Geneid'] == interval]
    PIC_invivo = PIC_invivo.append(df_temp2, ignore_index = True )
    
#We need to make bed files with the coordinates of the intervals, so we choose only the coordinates columns
Vibrio_in_vitro = Vibrio_invitro.iloc[:,0:5]
Vibrio_invivo = Vibrio_invivo.iloc[:,0:5]

PIC_invivo = PIC_invivo.iloc[:,0:5]


print(Vibrio_invitro)
print(Vibrio_invivo)
print(PIC_invivo)


#We don't differentially accessible peaks for PIC in-vitro, so there are none in the result files
Vibrio_invitro.to_csv('Vibrio_invitro.bed' , index = False)

Vibrio_invivo.to_csv('Vibrio_invivo.bed' , index = False)
PIC_invivo.to_csv('PIC_invivo.bed' , index = False)

###################################################################################################################################

