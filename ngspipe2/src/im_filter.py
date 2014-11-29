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
        print('Error: -i, -o and -p are required')
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
        outdir = args.o + '/' + 'imf'
    else:
        outdir = args.o + 'imf'
    ## Validate pedigree file??
    # create variant dataframe and pedigree dataframe
    df = pd.read_csv(args.i, header=0, sep='\t', dtype=str)
    ped_df = pd.read_csv(args.p, header=0, sep='\t', dtype=str) 
    # decide the family category
    print(outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
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
    if(imodel=='DeNovo'):
        print('DeNovo model: Filter variants by maf = ', config.denovo_maf)
        ngsdf = ngspipe_filter.filter_by_maf(ngsdf, config.denovo_maf)
    elif(imodel=='Recessive'):
        print('Recessive model: Filter variants by maf = ', \
              config.recessive_maf)
        ngsdf = ngspipe_filter.filter_by_maf(ngsdf, config.recessive_maf)
    elif(imodel=='Dominant'):
        print('Dominant model: Filter variants by maf = ', config.dominant_maf)
        ngsdf = ngspipe_filter.filter_by_maf(ngsdf, config.dominant_maf)
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
            filter_singleton(ngsdf, famdf, outdir, imodel)
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
    
def filter_singleton(ngsdf, famdf, outdir, imodel):
    if(imodel == 'DeNovo'):
        print('Filter singleton with DeNovo model ...')
        outdir = outdir + '/'+ imodel +'/singleton'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        sid = famdf.SampleID.tolist()[0]   
        outfile = outdir + '/' + sid+'.imfd'
             
        singletondf = ngsdf[config.im_header+[sid]]
        sngdf = singletondf.drop_duplicates(subset=config.vcf_basic_header)
        
        het_wm = generate_het_wm(config.allele_num)       
        sngdf = sngdf[sngdf[sid].isin(het_wm)]
        sngdf.to_csv(outfile, header=True, sep='\t', index=None)
        
    elif(imodel == 'Dominant'):
        print('Filter singleton with Dominant model ...')
    elif(imodel == 'Recessive'):
        print('Filter singleton with Recessive model ...')
        # find the homozygous or intrans based on vcf header
        outdir = outdir + '/'+ imodel +'/singleton'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        sid = famdf.SampleID.tolist()[0]   
        outfile = outdir + '/' + sid+'.imfd'
             
        singletondf = ngsdf[config.im_header+[sid]]
        sngdf = singletondf.drop_duplicates(subset=config.vcf_basic_header)
        
        homo_mm = generate_homo_mm(config.allele_num)
        sngdf = sngdf[sngdf[sid].isin(homo_mm)]
        sngdf.to_csv(outfile, header=True, sep='\t', index=None)
        
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
def generate_het_mm(allelenum): # return list
    het_mm = []
    for i in range(allelenum)[1:]:
        k=i+1
        for j in range(allelenum)[k:]:
            het_mm.append(str(i)+'/'+str(j))
            het_mm.append(str(i)+'|'+str(j))
    return het_mm
def generate_het_wm(allelenum):
    i = 0
    het_wm = []
    for j in range(allelenum)[1:]:
        het_wm.append(str(i)+'/'+str(j))
        het_wm.append(str(i)+'|'+str(j))
    return het_wm
def generate_homo_mm(allelenum):
    homo_mm=[]
    for i in range(allelenum)[1:]:
        homo_mm.append(str(i)+'/'+str(i))
        homo_mm.append(str(i)+'|'+str(i))
    return homo_mm
if __name__ == '__main__':
    main()