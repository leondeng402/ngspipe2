#! /usr/local/bin/python3.4
#####################################################################
#
#    ngspipe_filter.py: a program to filter the annotated ngs data
#    author:     Liyong Deng
#    copyright (c) 2014 by Liyong Deng in Dr. Chung's Lab
#
#####################################################################

import os
import argparse
import pandas as pd
import config


def main():   
    #####################################################################
    #    handle command line arguments
    #####################################################################
    parser = argparse.ArgumentParser(description='Annotate Next Generation ' \
                                     + 'Sequence file')
    parser.add_argument('-i', nargs='?', help='input file')
    parser.add_argument('-o', nargs='?', help='output directory')
    parser.add_argument('-f', nargs='?', help='filter file')
    args=parser.parse_args()
    if(args.i == None or args.o == None or args.f == None):
        print('Error: -i, -o and -f are required')
        exit()
    if(not os.path.isfile(args.i)):
        print('Error: File ' + args.i + ' does not exist.')
        exit()
    if(not args.i.endswith('.antd')):
        print('Error: file format should be .antd')
        exit()
    flt_df = prelimfilter(args.i, args.f)
    
    if(len(flt_df.index) != 0):
        infile_dir = os.path.dirname(args.i)
        infile_base = os.path.basename(args.i)
        outfile = infile_dir +'/' + infile_base.replace('.antd', '.fltd')
        flt_df.to_csv(outfile, header=True, index=False, sep='\t')
    else:
        print('No variatns left after filtering')
    
    
def prelimfilter(args_i, args_f): # return dataframe
    df = pd.read_csv(args_i, header=0, sep='\t', low_memory=False)
    header = list(df.columns.values)
    print('Start to filter variants ...')    
    maf = 0
    vflist=[]
    eflist=[]
    pathogenic_params ={}    
    print('Parse filter file ...')
    with open(args_f, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if(line == '' or line.startswith('#')):
                continue
            formula = line.split('&')        
            if(len(formula) == 1):
                words=formula[0].split('=')
                if(words[0] == 'maf'):
                    maf=float(words[1])
                elif(words[0] == 'VariantFunction'):
                    vflist=words[1].split(',')
                    vflist = [x.strip() for x in vflist]
                elif(words[0] == 'ExonicFunction'):
                    eflist=words[1].split(',')
                    eflist = [x.strip() for x in eflist]
                else:
                    print(words[0])
                    print('Error: not correct grammar for a filter file')
                    return pd.DataFrame(columns=header)
            elif(len(formula)== 2):
                words0 = formula[0].split('=')
                words1 = formula[1].split('=')
                words0[0]=words0[0].strip()
                words0[1]=words0[1].strip()        
                if(words0[0]=='ExonicFunction' and words0[1] == \
                   'nonsynonymous SNV'):
                    words1[0]=words1[0].strip()
                    words1[1]=words1[1].strip()
                    pathogenic_params[words1[0]]=words1[1]
                
                else: 
                    print(formula, words0[0], words0[1], '2')
                    print('Error: not correct grammar for a filter file')
                    return pd.DataFrame(columns=header)
            else:
                print(formula, '>2')
                print('Error: not correct grammar for a filter file')
                return pd.DataFrame(columns=header)
    print('Start to filter by maf ...')    
    df = filter_by_maf(df, maf)
    
    print('Start to filter by VariantFunction ...')
    df = filter_by_vf(df, vflist)
    
    print('Start to filter by ExonicFunction ...')
    df = filter_by_ef(df, eflist)
    
    print('Start to filter nonsynoymous snv by pathogenicity ...')
    ns_snv_df=df[df['ExonicFunction']=='nonsynonymous SNV']
    rest_df=df[df['ExonicFunction']!='nonsynonymous SNV']    
    ns_snv_filtered_df = filter_by_pathogenicity(ns_snv_df, pathogenic_params)    
    df = rest_df.append(ns_snv_filtered_df, ignore_index=True)
        
    return df

def filter_by_ef(dataframe, eflist): # return DataFrame
    header=list(dataframe.columns.values)
    if(len(eflist) <= 0):
        return pd.DataFrame(columns=header)
    temp_df = dataframe[dataframe.ExonicFunction.isin(eflist)]     
    return temp_df
   
def filter_by_maf(dataframe, maf): # return DataFrame
    rsvd_rows = []
    dump_rows = []
    header = list(dataframe.columns.values)
    af_header= config.allele_frequency_header
    dataframe[af_header] = dataframe[af_header].astype(float)
    for index, row in dataframe.iterrows():
        af_list = row[config.allele_frequency_header].tolist()        
        max_af = max(af_list)        
        if(max_af <= maf):
            rsvd_rows.append(row)
        else:
            dump_rows.append(row)
    rsvd_df = pd.DataFrame(rsvd_rows, columns=header)
    dump_df = pd.DataFrame(dump_rows, columns=header)
    return rsvd_df

def filter_by_pathogenicity(dataframe, pathogenic_params):# return DataFrame
    header=list(dataframe.columns.values)
    df = pd.DataFrame(columns=header)
    for key, value in pathogenic_params.items():           
        key=key.strip()
        value=value.strip()
        if('CADD' in key):                          
            temp_df = dataframe[dataframe[key]>=float(value)]  
            dataframe = dataframe[dataframe[key]<float(value)]             
            df=df.append(temp_df, ignore_index=True)            
        else:
            # sift and polyphen2 
            temp_df = dataframe[dataframe[key]==value] 
            dataframe= dataframe[dataframe[key]!=value] 
            df=df.append(temp_df, ignore_index=True)  
    df[header]=df[header].astype(str)
    return df      
    
def filter_by_vf(dataframe, vflist): # return DataFrame
    
    header=list(dataframe.columns.values)
    if(len(vflist) <= 0):
        return pd.DataFrame(columns=header)
    temp_df = dataframe[dataframe.VariantFunction.isin(vflist)]    
    return temp_df


if __name__ == '__main__':
    main()