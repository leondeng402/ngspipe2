#! /usr/local/bin/python3.4
#####################################################################
#
#    ino.py: a module contains methods extracting information from 
#            input file and output information in a certain format
#    copyright (c) 2014 by Liyong Deng in Dr. Chung's Lab
#
#####################################################################
import pandas as pd
import config
# All functions organized in alphabetic order

def get_genotype(inputfile, ids): #return a DataFrame
    results=[]
    with open(inputfile, 'r')as f:        
        for line in f.readlines():
            line = line.strip()    
            words=line.split('\t')
            variant=words[:5]
            result=words[8:13]+words[17:]
            row=[]
            for words in result:                    
                word=words.split(':')
                row.append(word[0])
            results.append(variant+row)
    header = config.basic_header + config.vcf_basic_header + ids 
    genotype_df = pd.DataFrame(results,  columns=header)
    #genotype_df = pd.DataFrame(results, index=None, columns=header, \
    #                           dtype='str')
    return genotype_df
    
def get_sampleid(inputfile): #return a list
    sampleids=None
    with open(inputfile, 'r')as f:        
        for line in f.readlines():
            line = line.strip()    
            if(line.startswith('#CHROM')):            
                line = line.replace('#', '')
                words=line.split('\t')                
                sampleids=words[9:]
                break
    return sampleids
