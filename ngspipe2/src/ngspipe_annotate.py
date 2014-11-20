#! /usr/local/bin/python3.4
#####################################################################
#
#    ngspipe_annotate.py: a program to annotate next generation sequencing 
#                data at variant,gene and region levels
#    author:     Liyong Deng
#    copyright (c) 2014 by Liyong Deng in Dr. Chung's Lab
#
#####################################################################

import os, sys, time, subprocess, shlex
import argparse

import pandas as pd
import config, ino
#import config, gene, variant, region, operation

def main(args):   
    # Validate the arguments
    if(args.i == None):
        print('Error: -i is required')
        exit()
    if(not os.path.isfile(args.i)):
        print('Error: File ' + args.i + ' does not exist.')
        exit()
    if(args.o!= None and not os.path.isdir(args.o)):
        print('Error: Directory ' + args.o + ' does not exist.')
        exit()
    infile = args.i
    # time the process
    if(args.t):
        start_time = time.time()
    # Two types of input file: txt and vcf
    print('\n**************************************************************')
    print('*    Parse the input file and Prepare the anvinput file      *')
    print('**************************************************************')
    vcf_flag = False
    fieldnum = 0
    
    if(infile.endswith('.vcf')):
        print('This is vcf file')
        vcf_flag = True
        infile_dir = os.path.dirname(infile)+'/'
        infile_base = os.path.basename(infile)
        # extract sample ids in #CHROM line
        sampleids = ino.get_sampleid(infile)    
        # prepare anvinput and retrieve genotype
        outfile = infile_dir + infile_base.replace('.vcf', '.temp')
        md_str = 'convert2annovar.pl  -format vcf4 -allsample -withfreq  '\
                 + '-include ' + infile + ' -outfile ' + outfile
        cmd = shlex.split(md_str)
        p = subprocess.Popen(cmd)
        p.wait()
        
        ## retrieve genotypes
        geno_df = ino.get_genotype(outfile, sampleids)
        os.remove(outfile)
        ## create the avinput file 
        avinput = outfile.replace('.temp', '.txt')
        avinput_header = config.basic_header
        geno_df.to_csv(avinput, header=False, sep='\t', index=False, \
                   columns=avinput_header)  
        fieldnum = len(avinput_header) 
    elif(infile.endswith('.txt')):
        print('This is txt file')
        # decide fieldnum in txt file
        for line in open(infile):
            words = line.split('\t')
            fieldnum = len(words)
            break
        
        if(fieldnum not in range(5, 7)):    
            print('Error: columns number ' + str(fieldnum) + \
                  '  is out of range[5,6]')
            exit()
        avinput = infile
    else:
        print('Error: input file has a wrong format')
    
    # time the process
    if(args.t):
        end_time = time.time()
        print('Time consumed: ', \
              time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
        start_time = time.time()
        
    #####################################################################
    #  table_annovar queries in one
    #####################################################################
    print('\n**************************************************************')
    print('*    multiple database annotations                           *')
    print('**************************************************************')
    md_str= 'table_annovar.pl ' + avinput + ' ' + config.dir_humandb \
            + ' -buildver ' + config.buildver \
            + ' -protocol ' + ','.join(config.protocols) \
            + ' -operation ' + ','.join(config.operations)\
            + ' -nastring ' + config.nastring \
            + ' -csvout'
   
    if(args.r):
        md_str=md_str + ' -remove'
        #print(md_str)
    cmd = shlex.split(md_str)
    subprocess.call(cmd)
    result = avinput+'.'+config.buildver + config.tblanno_suffix
    header=config.tblanno_header
    df = pd.read_csv(result, header=None, names=header, sep=',', \
                 dtype='str', low_memory=False, skiprows=1)
    df = df.fillna(config.nastring)
    print(df.iloc[0:5])
    # time the process
    if(args.t):
        end_time = time.time()
        print('Time consumed: ', \
              time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
        start_time = time.time()
    return 0


if __name__ == '__main__':
    #####################################################################
    #    handle command line arguments
    #####################################################################
    parser = argparse.ArgumentParser(description='Annotate Next Generation ' \
                                     + 'Sequence file')
    parser.add_argument('-e', action='store_true', help='expanding gene ' \
                        + 'entries')
    parser.add_argument('-i', nargs='?', help='input file')
    parser.add_argument('-m', action='store_true', help='merge gene info ' \
                    + 'from refgene, knowngene, and ensgene')
    parser.add_argument('-o', nargs='?', help='output directory')
    parser.add_argument('-r', action='store_true', help='delete temporary '\
                        +'files')
    parser.add_argument('-t', action='store_true', help='timing each '\
                        +'operation')
    parser.add_argument('-p', nargs='?', help='pedigree file')

    arguments=parser.parse_args()
    
    main(arguments)
    
    