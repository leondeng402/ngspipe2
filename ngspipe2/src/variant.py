#! /usr/local/bin/python3.4
#####################################################################
#
#    variant.py: a module contains methods to do variant level 
#                annotation and etc
#    author:  Liyong Deng
#    copyright (c) 2014 by Liyong Deng in Dr. Chung's Lab
#
#####################################################################
import os
import subprocess, shlex
import pandas as pd
import config
def fetch_cadd(inputfile, fieldnum, vflag, rflag): # return Datsasame    
    ### fetch popfreq_allele_freq
    print('\n**************************************************************') 
    print('*    Fetch cadd pathogenicity scores                         *')
    print('**************************************************************')
    cmd_str= 'annotate_variation.pl  -filter -otherinfo -out ' + inputfile   \
             + ' -dbtype ' + config.cadd + ' -build ' + config.buildver + ' ' \
             + inputfile  + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.cadd_drp_suffix
    of = open(result+'.temp', 'w')
    with open(result, 'r') as f:
        data = f.read().replace('\t', ',')
        of.write(data)    
    of.close()
    os.remove(result)
    os.rename(result+'.temp', result)
    if(vflag or fieldnum == 5):
        header = config.cadd_header
    else:
        header = config.cadd_header[:3] + config.txtinput_header
    cadd_df = pd.read_csv(result, header=None, names=header, \
                             sep=',', usecols=header[1:],\
                             low_memory=False)
    if 'Comments' in cadd_df:
        cadd_df= cadd_df.drop('Comments', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.cadd_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.cadd_flt_suffix)
    if (len(cadd_df.index) == 0):
        cadd_df = pd.DataFrame(columns=config.cadd_use_header)
    cadd_df[config.basic_header] = cadd_df[config.basic_header].astype(str)
    return cadd_df

def fetch_hgmd(): # return DataFrame
    print('\n**************************************************************')
    print('*    Start to fetch hgmd                                     *')
    print('**************************************************************')
    hgmd_df = pd.read_csv(config.hgmd_file, header= None, 
                          names=config.hgmd_header, \
                          sep='!', index_col=False, 
                          usecols=config.hgmd_use_header, \
                          low_memory=False)
    hgmd_df[['Chr', 'Start', 'End']] = hgmd_df[['Chr', 'Start', 'End']].\
                                       astype(str)
    return hgmd_df

def fetch_g1k(inputfile, race, fieldnum, vflag, rflag): # return DataFrame    
    ### fetch g1k_all_allele_freq
    print('\nFetch g1k_' + race + ' allele frequency')
    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile  \
             + ' -dbtype ' + config.g1k + '_' + race \
             +' -build ' + config.buildver + ' '\
             + inputfile  + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+'_' + race.upper() \
             + config.g1k_drp_suffix
    if(vflag or fieldnum==5):
        header = config.g1k_header
    else:
        header = config.g1k_header[:2] + config.txtinput_header
    g1k_df = pd.read_csv(result, header=None, names=header, \
                             sep='\t', usecols=header[1:],\
                             low_memory=False)
    g1k_df = g1k_df.rename(columns={'Value': config.g1k.upper() + '_' \
                                    + race.upper()})
    if(rflag):
        os.remove(inputfile+'.'+config.buildver + '_' + race.upper() \
                  + config.g1k_drp_suffix)
        os.remove(inputfile+'.'+config.buildver + '_' + race.upper() \
                  +config.g1k_flt_suffix)
    if (len(g1k_df.index) == 0):
        g1k_df = pd.DataFrame(columns=config.g1k_header[1:])
        g1k_df = g1k_df.rename(columns={'Value': config.g1k.upper() + '_' \
                                        + race.upper()})

    if 'Comments' in g1k_df:
        g1k_df= g1k_df.drop('Comments', axis=1)
    g1k_df[config.basic_header] = g1k_df[config.basic_header].astype(str)

    return g1k_df

def fetch_esp(inputfile, race, fieldnum, vflag, rflag): # return DataFrame    
    ### fetch g1k_all_allele_freq
    print('\nFetch esp_' + race + ' allele frequency')
    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile  \
             + ' -dbtype ' + config.esp + '_' + race \
             +' -build ' + config.buildver + ' '\
             + inputfile  + ' ' + config.dir_humandb
    print(cmd_str)
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+'_' + config.esp  + '_' \
            + race + config.esp_drp_suffix
    print(result)
    if(vflag or fieldnum==5):
        header = config.esp_header
    else:
        header = config.esp_header[:2] + config.txtinput_header
  
    esp_df = pd.read_csv(result, header=None, names=header, \
                             sep='\t', usecols=header[1:],\
                             low_memory=False)
    esp_df = esp_df.rename(columns={'Value': config.esp.upper() + '_' \
                                    + race.upper()})
    if(rflag):
        os.remove(inputfile+'.'+config.buildver + '_' + config.esp + '_' \
                  + race + config.esp_drp_suffix)
        os.remove(inputfile+'.'+config.buildver + '_' + config.esp + '_' \
                  + race + config.esp_flt_suffix)
    if (len(esp_df.index) == 0):
        esp_df = pd.DataFrame(columns=config.esp_header[1:])
        esp_df = esp_df.rename(columns={'Value': config.esp + '_' \
                                        + race})

    if 'Comments' in esp_df:
        esp_df= esp_df.drop('Comments', axis=1)
    esp_df[config.basic_header] = esp_df[config.basic_header].astype(str)
    
    return esp_df