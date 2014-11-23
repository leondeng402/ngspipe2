#! /usr/local/bin/python3.4
#####################################################################
#
#    gene.py:  the file contains all the gene level functions#             
#    author:  Liyong Deng
#    copyright (c) 2014 by Liyong Deng
#
#####################################################################
import os
import subprocess, shlex

import pandas as pd

import config

#####################################################################
# expand_gene_entries
#####################################################################
def expand_gene_entries(dataframe):
    print('\n**************************************************************')
    print('*    Expanding gene entries                                  *')
    print('**************************************************************')   
    headers = list(dataframe.columns.values)    
    temp_list = []
    j=0
    for index, row in dataframe.iterrows():
        #iterate the dataframe to find the row with multiple gene entries
        gene_piece = row['Gene'].strip().split(',')
        aa_change_piece = row['AAChange'].strip().split(',')
        for gene_item in gene_piece:
            for aa_item in aa_change_piece:
                #print('gene: ', gene_piece)
                temp_row=pd.Series(row.values.copy(), index=row.index, \
                                   name=row.name)
                temp_row['Gene'] = gene_item
                temp_row['AAChange']= aa_item
                temp_list.append(temp_row)           
    temp_df = pd.DataFrame(temp_list, columns=headers)            
    return temp_df

def merge_gene_entries(dataframe): # return a DataFrame
    #  refgene query, save known to gene_df and unknown to unknown_df
    #headers = config.merged_custom_header
    headers = list(dataframe.columns.values)
    #print(headers)
    # remove the knowngene and ensgene headers
    for item in config.knowngene_header + config.ensgene_header:
        #print(item)
        headers.remove(item)
    if(  config.refgene_header[0] not in headers):
        print('Error: refGene annotation is not in')
        exit()
    print('\n**************************************************************')
    print('*    merging gene entries from refgene, knowngene, & ensgene *')
    print('**************************************************************')      
    temp_list = []
    j=0
    for index, row in dataframe.iterrows():
        temp_list.append(row[headers])
        # check knowngene gene entries
        temp_row=pd.Series(row.values.copy(), index=row.index, \
                                   name=row.name)
        temp_row = temp_row[headers]        
        if(row['VariantFunction'] != row['VariantFunction.k'] \
           or row['ExonicFunction'] != row['ExonicFunction.k']) :
            #print('refGene & knownGene info is not the same')
            #print(temp_row['Chr'], temp_row['Start'], temp_row['End'], \
            #      temp_row['Ref'], temp_row['Alt'])
            temp_row['VariantFunction'] = row['VariantFunction.k']
            temp_row['Gene'] = row['Gene.k']
            temp_row['GeneDetail'] = row['GeneDetail.k']
            temp_row['ExonicFunction'] = row['ExonicFunction.k']
            temp_row['AAChange'] = row['AAChange.k']
            temp_list.append(temp_row)
        # check ensngene gene entries
        temp_row=pd.Series(row.values.copy(), index=row.index, \
                                   name=row.name)
        temp_row = temp_row[headers]       
        if(row['VariantFunction'] != row['VariantFunction.e'] \
           or row['ExonicFunction'] != row['ExonicFunction.e']):
            #print('refGene & ensGene info is not the same')
            #print(temp_row['Chr'], temp_row['Start'], temp_row['End'], \
            #      temp_row['Ref'], temp_row['Alt'])
            temp_row['VariantFunction'] = row['VariantFunction.e']
            temp_row['Gene'] = row['Gene.e']
            temp_row['GeneDetail'] = row['GeneDetail.e']
            temp_row['ExonicFunction'] = row['ExonicFunction.e']
            temp_row['AAChange'] = row['AAChange.e']
            temp_list.append(temp_row)
    temp_df = pd.DataFrame(temp_list, columns=headers)
            
    return temp_df