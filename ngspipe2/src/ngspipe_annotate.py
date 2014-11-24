#! /usr/local/bin/python3.4
#####################################################################
#
#    ngspipe_annotate.py: a program to annotate next generation sequencing 
#                data at variant,gene and region levels
#    author:     Liyong Deng
#    copyright (c) 2014 by Liyong Deng in Dr. Chung's Lab
#
#####################################################################

import os, time, subprocess, shlex
import argparse
import pandas as pd
import config, ino, gene, variant

    
def main():   
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
    args=parser.parse_args()
    # Validate the arguments
    if(args.i == None or args.o == None):
        print('Error: -i and -o are required')
        exit()
    if(not os.path.isfile(args.i)):
        print('Error: File ' + args.i + ' does not exist.')
        exit()
    if(args.o!= None and not os.path.isdir(args.o)):
        print('Error: Directory ' + args.o + ' does not exist.')
        exit()
    annotate(args.e, args.i, args.m, args.o, args.p, args.r, args.t)
        
def annotate(args_e, args_i, args_m, args_o, args_p, args_r, args_t):
    infile = args_i
    # time the process
    if(args_t):
        start_time = time.time()
    # Two types of input file: txt and vcf
    print('\n**************************************************************')
    print('*    Parse the input file and Prepare the anvinput file      *')
    print('**************************************************************')
    vcf_flag = False
    fieldnum = 0
    infile_dir = os.path.dirname(infile)
    infile_base = os.path.basename(infile)
    
    if(infile.endswith('.vcf')):
        print('This is vcf file')
        vcf_flag = True
        
        # extract sample ids in #CHROM line
        sampleids = ino.get_sampleid(infile)    
        # prepare anvinput and retrieve genotype
        outfile = infile_dir +'/'+ infile_base.replace('.vcf', '.temp')
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
        return None
    
    # time the process
    if(args_t):
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
   
    if(args_r):
        md_str=md_str + ' -remove'
    cmd = shlex.split(md_str)
    subprocess.call(cmd)
    result = avinput+'.'+config.buildver + config.tblanno_suffix
    header=config.tblanno_header
    df = pd.read_csv(result, header=None, names=header, sep=',', \
                 dtype='str', low_memory=False, skiprows=1)
    df = df.fillna(config.nastring)
    
    # time the process
    if(args_t):
        end_time = time.time()
        print('Time consumed: ', \
              time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
        start_time = time.time()
        
    #####################################################################
    #  fetch 1000 genome project all
    #####################################################################
    for race in config.g1k_races:
        g1k_df = variant.fetch_g1k(avinput, race, fieldnum, vcf_flag, args_r)
        
        #print(g1k_df)
        df = df.merge(g1k_df, how='left', on=config.basic_header)
    df = df.fillna(0)

    if(args_t):
        end_time = time.time()
        print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
        start_time = time.time()
    #####################################################################
    #  fetch esp project allele frequency
    #####################################################################
    for race in config.esp_races:
        esp_df = variant.fetch_esp(avinput, race, fieldnum, vcf_flag, args_r)
        df = df.merge(esp_df, how='left', on=config.basic_header)
    df = df.fillna(0)

    if(args_t):
        end_time = time.time()
        print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
        start_time = time.time()   
    #####################################################################
    #  fetch cadd score
    #####################################################################
    cadd_df = variant.fetch_cadd(avinput, fieldnum, vcf_flag, args_r)
    
    #print(cadd_df)
    df = pd.merge(df, cadd_df, how='left', on=config.basic_header)
    df = df.fillna(0)
    if(args_r):
        os.remove(avinput+'.log')
    
    if(args_t):
        end_time = time.time()
        print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
        start_time = time.time()
    #####################################################################
    #  fetch hgmd info
    #####################################################################
    hgmd_df = variant.fetch_hgmd()    
    df = df.merge(hgmd_df, how='left', on=['Chr', 'Start', 'End'])
    df = df.fillna('.')
    if(args_t):
        end_time = time.time()
        print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
        start_time = time.time()
        #print(df.columns.values)
    if(vcf_flag):
        #####################################################################
        #  add genotypes to annotated variants
        #####################################################################
        print('\n**************************************************************')
        print('*    add genotype info to annotated variants                 *')
        print('**************************************************************')
        df = df.merge(geno_df, how='left', on=config.basic_header)
        if(args_t):
            end_time = time.time()
            print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    #####################################################################
    #  merge gene entries from different sources
    #####################################################################
    print(args_o, infile_base)
    if(args_o.endswith('/')):
        outfile = args_o + infile_base + '.antd'
    else:
        outfile = args_o + '/' + infile_base + '.antd'
    if(args_m):
        # merging gene entries
        df = gene.merge_gene_entries(df) 
        outfile = outfile.replace('.antd', '.mrgd.antd')      
        df.to_csv(outfile, header=True, sep='\t', index=False)
        if(args_t):
            end_time = time.time()
            print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    
    #####################################################################
    #  expand multiple entries in gene info
    #####################################################################
    if(args_e):
        # expanding gene entries
        df = gene.expand_gene_entries(df) 
        outfile = outfile.replace('.antd', '.xpnd.antd') 
        df.to_csv(outfile, header=True, sep='\t', index=False)
        if(args_t):
            end_time = time.time()
            print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    
    if(not args_m and not args_e):
        df.to_csv(outfile, header=True, sep='\t', index=False,)   
    return df


if __name__ == '__main__':
    main()

    
    