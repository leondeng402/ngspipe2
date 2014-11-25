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
    parser.add_argument('-m', nargs='?', choices=['DeNovo', 'Dominant', \
                                                  'Recessive'],\
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
    if(not (args.i.endswith('.fltd') or  args.i.endswith('.rnkd'))):
        print('Error: file format should be either .fltd or .rnkd')
        exit()
    if(args.o!= None and not os.path.isdir(args.o)):
        print('Error: Directory ' + args.o + ' does not exist.')
        exit()
    if(not args.o.endswith('/')):
        outdir = args.o + '/'
    else:
        outdir = args.o
    ## Validate pedigree file??
    # create variant dataframe and pedigree dataframe
    df = pd.read_csv(args.i, header=0, sep='\t', low_memory=False)
    ped_df = pd.read_csv(args.p, header=0, sep='\t', low_memory=False) 
    # decide the family category
    filter_by_family(df, ped_df, outdir, args.m)
    
def filter_by_dominant(ngsdf, peddf, outdir): 
    pass
    # filter by maf
    
    # proband are hets while both affected parents are hets 
    
    
def filter_by_recessive(ngsdf, peddf, outdir): 
    pass   
    # filter by maf
    
    # 1. filter by homo alt alleles in proband and unaffected parents are het 
    
    # 2. filter by at least two hets in the same gene and parents carry each
    # alt allele respectively
    
           
def filter_by_family(ngsdf, peddf, outdir, imodel):    
    # filter by maf
    ngsdf = ngspipe_filter.filter_by_maf(ngsdf, config.denovo_maf)
    
    # decide singleton,duos, sngtnfam, trios, duofam, triofam, other
    probandids = peddf[peddf.Person=='Proband'].SampleID.tolist()
    for pid in probandids:
        familyid = peddf[peddf.SampleID == pid].FamilyID.tolist()
        famdf = peddf[peddf.FamilyID.isin(familyid)]
        
        #print(familyid, famdf)
        family_size = len(famdf.index)
        print('\n'+familyid[0]+': ', family_size)
        if(family_size == 1):
            print('This is a singleton')
            filter_singleton(famdf, peddf, outdir, imodel)
        elif(family_size == 2):
            personlist = famdf.Person.tolist()
            if('Mother' in personlist or 'Fatehr' in personlist):
                print('This is a duos')
                filter_duos(famdf, peddf, outdir, imodel)
            else:
                print('This is sngtnfam')
                filter_family_with_singleton(famdf, peddf, outdir, imodel)
            
        elif(family_size == 3):
            personlist = famdf.Person.tolist()
            #print(personlist)
            if(('Mother' in personlist) and ('Father' in personlist)):
                print('This is a trios')
                filter_trios(famdf, peddf, outdir, imodel)
            elif(('Mother' in personlist) != ('Father' in personlist)):
                print('This a duofam')
                filter_family_with_duos(famdf, peddf, outdir, imodel)
            else:
                print('This is a sngtnfam')
                filter_family_with_singleton(famdf, peddf, outdir, imodel)
                
        elif(family_size >=4):
            personlist = famdf.Person.tolist()
            if(('Mother' in personlist) and ('Father' in personlist)):
                print('This is a triofam')
                filter_family_with_trios(famdf, peddf, outdir, imodel)
            elif(('Mother' in personlist) != ('Father' in personlist)):
                print('This a duofam')
                filter_family_with_duos(famdf, peddf, outdir, imodel)
            else:
                print('This is a sngtnfam')
                filter_family_with_singleton(famdf, peddf, outdir, imodel)
    # filter by het in proband
    
    # filter by that alt allele in probnad doesn't exist in parents
    
def filter_singleton(ngsdf, peddf, outdir, imodel):
    if(imodel == 'DeNovo'):
        print('Filter singleton with DeNovo model ...')
    elif(imodel == 'Dominant'):
        print('Filter singleton with Dominant model ...')
    elif(imodel == 'Recessive'):
        print('Filter singleton with Recessive model ...')
def filter_family_with_singleton(ngsdf, peddf, outdir, imodel):
    if(imodel == 'DeNovo'):
        print('Filter family with singleton with DeNovo model ...')
    elif(imodel == 'Dominant'):
        print('Filter family with singleton with Dominant model ...')
    elif(imodel == 'Recessive'):
        print('Filter family with singleton with Recessive model ...')
def filter_duos(ngsdf, peddf, outdir, imodel):
    if(imodel == 'DeNovo'):
        print('Filter duos with DeNovo model ...')
    elif(imodel == 'Dominant'):
        print('Filter duos with Dominant model ...')
    elif(imodel == 'Recessive'):
        print('Filter duos with Recessive model ...')
def filter_family_with_duos(ngsdf, peddf, outdir, imodel):
    if(imodel == 'DeNovo'):
        print('Filter family with duos with DeNovo model ...')
    elif(imodel == 'Dominant'):
        print('Filter family with duos with Dominant model ...')
    elif(imodel == 'Recessive'):
        print('Filter family with duos with Recessive model ...') 
def filter_family_with_trios(ngsdf, peddf, outdir, imodel):
    if(imodel == 'DeNovo'):
        print('Filter family with trios with DeNovo model ...')
    elif(imodel == 'Dominant'):
        print('Filter family with trios with Dominant model ...')
    elif(imodel == 'Recessive'):
        print('Filter family with trios with Recessive model ...')      
def filter_trios(ngsdf, peddf, outdir, imodel):
    if(imodel == 'DeNovo'):
        print('Filter trios with DeNovo model ...')
    elif(imodel == 'Dominant'):
        print('Filter trios with Dominant model ...')
    elif(imodel == 'Recessive'):
        print('Filter trios with Recessive model ...')   
if __name__ == '__main__':
    main()