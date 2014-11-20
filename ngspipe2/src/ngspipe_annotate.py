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

#import config, gene, variant, region, operation

def main(args):
   
    
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
    # Two types of input file: txt and vcf
    vcf_flag = False
    fieldnum = 0
    if(infile.endswith('.vcf')):
        print('This is vcf file')
        vcf_flag = True
        infile_dir = os.path.dirname(infile)+'/'
        infile_base = os.path.basename(infile)
        # create anvinput and retrieve genotype
        outfile= infile_dir + infile_base.replace('.vcf', '.temp')
        md_str = 'convert2annovar.pl  -format vcf4 -allsample -withfreq  '\
                 + '-include ' + infile + ' -outfile ' + outfile
        cmd = shlex.split(md_str)
        p = subprocess.Popen(cmd)
        p.wait()
    elif(infile.endswith('.txt')):
        print('This is txt file')
        # decide fieldnum in txt file
        for line in open(infile):
            words = line.split('\t')
            fieldnum = len(words)
            break
    
        if(fieldnum not in range(5, 7)):    
            print('Error: columns number ' + str(fieldnum) + '  is out of range[5,6]')
            exit()
    else:
        print('Error: input file has a wrong format')
    
    
    
    
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
    print(arguments)
    arguments.i ='b'
    main(arguments)
    
    