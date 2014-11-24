#! /usr/local/bin/python3.4
#####################################################################
#
#    im_filter.py: a program to filter the annotated ngs data using
#                inheritance models
#    author:     Liyong Deng
#    copyright (c) 2014 by Liyong Deng in Dr. Chung's Lab
#
#####################################################################

import os
import argparse
import pandas as pd
import config, ngspipe_filter

def main():
    #####################################################################
    #    handle command line arguments
    #####################################################################
    parser = argparse.ArgumentParser(description='Filter annotated next ' \
                                     + 'generation sequence file using ' \
                                     +'inheritance models')
    parser.add_argument('-i', nargs='?', help='input file')
    parser.add_argument('-m', nargs='?', choices=['DeNovo', 'AD', 'AR', 'XR'],\
                    default='DeNovo',  help='inheritance models to be used')
    parser.add_argument('-o', nargs='?', help='output directory')
    parser.add_argument('-p', nargs='?', help='pedigree file')
    args=parser.parse_args()
    if(args.i == None or args.o == None or args.p == None):
        print('Error: -i, -o and -f are required')
        exit()
    if(not os.path.isfile(args.i)):
        print('Error: File ' + args.i + ' does not exist.')
        exit()
    if(not args.i.endswith('.fltd') or not args.i.endswith('.rnkd')):
        print('Error: file format should be either .antd or .rnkd')
        exit()
    if(args.o!= None and not os.path.isdir(args.o)):
        print('Error: Directory ' + args.o + ' does not exist.')
        exit()
          
    ## Validate pedigree file??
    # create variant dataframe and pedigree dataframe
    df = pd.read_csv(args.i, header=0, sep='\t', low_memory=False)
    ped_df = pd.read_csv(args.p, header=0, sep='\t', low_memory=False) 
    # filter_by_model
    if(args.m == 'DeNovo'):
        filter_by_denovo(df, ped_df, args.o)
    elif(args.m == 'AD'):
        pass
    elif(args.m == 'AR'):
        pass
    elif(args.m == 'XR'):
        pass
    else:
        print('Error: wrong inheritance model, model ' + args.m + ' is not ' \
              +'in [DeNovo, AD, AR, XR]')
        exit()
    
def filter_by_ad(ngsdf, peddf, outdir): 
    pass
    # filter by maf
    
    # proband are hets while both affected parents are hets 
    
    
def filter_by_ar(ngsdf, peddf, outdir): 
    pass   
    # filter by maf
    
    # 1. filter by homo alt alleles in proband and unaffected parents are het 
    
    # 2. filter by at least two hets in the same gene and parents carry each
    # alt allele respectively
    
           
def filter_by_denovo(ngsdf, peddf, outdir):    
    # filter by maf
    ngsdf = ngspipe_filter.filter_by_maf(ngsdf, config.denovo_maf)
    
    # singleton/proband
    # duo
    # trios
    # family
    
    # filter by het in proband
    
    # filter by that alt allele in probnad doesn't exist in parents
    
def filter_by_xr(ngsdf, peddf, outdir):
    pass
    
    # filter by maf
    
    # 1. filter by homo alt alleles in proband and unaffected parents are het 
    
    # 2. filter by at least two hets in the same gene and parents carry each
    # alt allele respectively

def filter_by_het(dataframe, sampleid): 
    pass
if __name__ == '__main__':
    main()