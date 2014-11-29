#! /usr/local/bin/python3.4
#####################################################################
#
#    config.py: the configuration file for database files and their
#               headers, and file names for the intermediate files
#    author:  Liyong Deng
#    copyright (c) 2014 by Liyong Deng in Dr. Chung's Lab
#
#####################################################################
buildver = 'hg19'
dir_humandb = '/home/leon/biobin/annovar20141017/annovar/humandb'
basic_header = ['Chr', 'Start', 'End', 'Ref', 'Alt']
vcf_basic_header = ['CHROM', 'POS', 'ID', 'REF', 'ALT']
vcf_other_header = ['QUAL', 'FILTER', 'INFO', 'FORMAT']
txtinput_header = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Comments'] 
avinput_header = basic_header +['A', 'B','C'] + vcf_basic_header \
                 + vcf_other_header
                 
#####################################################################
#  region-based: -regionanno
#####################################################################
cytoband = 'cytoBand'
superdup = 'genomicSuperDups'
gwas = 'gwasCatalog'

#####################################################################
#  gene level:'--geneanno -dbtype refGene'
#####################################################################
refgene = 'refGene'
knowngene='knownGene'
ensgene='ensGene'
refgene_header = ['VariantFunction', 'Gene', 'GeneDetail', 'ExonicFunction', \
                  'AAChange']
knowngene_header = ['VariantFunction.k', 'Gene.k', 'GeneDetail.k', \
                    'ExonicFunction.k', 'AAChange.k']
ensgene_header = ['VariantFunction.e', 'Gene.e', 'GeneDetail.e', \
                  'ExonicFunction.e', 'AAChange.e']
####################################################################
#  variant leve: filter-based
####################################################################
snp='snp138'
cg='cg69'
#popfreq='popfreq_all'
ljb='ljb23_all'
clinvar='clinvar_20140929'
cosmic='cosmic70'
nci='nci60'
# esp allele frequency
esp='esp6500si'
esp_drp_suffix='_dropped'
esp_flt_suffix='_filtered'
esp_races=['all', 'aa', 'ea']
esp_header=['Label', 'Value', 'Chr', 'Start', 'End', 'Ref', 'Alt']
esp_use_header = [esp.upper() + '_' + esp_races[0].upper(), \
                  esp.upper() + '_' + esp_races[0].upper(), \
                  esp.upper() + '_' + esp_races[0].upper()]
# g1k allele frequency
g1k='1000g2014sep'
g1k_drp_suffix='.sites.2014_09_dropped'
g1k_flt_suffix='.sites.2014_09_filtered'
g1k_races=['all', 'afr', 'amr', 'eas', 'eur', 'sas']
g1k_header=['Label', 'Value', 'Chr', 'Start', 'End', 'Ref', 'Alt']
g1k_use_header = [g1k.upper() + '_' + g1k_races[0].upper(), \
                  g1k.upper() + '_' + g1k_races[1].upper(), \
                  g1k.upper() + '_' + g1k_races[2].upper(), \
                  g1k.upper() + '_' + g1k_races[3].upper(), \
                  g1k.upper() + '_' + g1k_races[4].upper(), \
                  g1k.upper() + '_' + g1k_races[5].upper()]

########################
# exac allele frequency
########################


# cadd raw score and phred score
cadd='caddgt10'
cadd_header = ['Label','CADD_raw', 'CADD_phred', 'Chr', 'Start', 'End', \
               'Ref', 'Alt']
cadd_use_header = cadd_header[1:]
cadd_drp_suffix='_'+cadd+'_dropped'
cadd_flt_suffix='_'+cadd+'_filtered'
# hgmd database configuration
hgmd_file= '/home/leon/lab/hgmd/hgmd_pro.allmut.txt'
hgmd_header = ['HGMD_Disease', 'HGMD_Gene', 'Chrom', 'HGMD_Genename', 'gdbid', 
               'Omimid', 'Amino', 'Deletion', 'Insertion', 'HGMD_Codon', 
               'CodonAff', 'Descr','Hgvs', 'HGMD_hgvsAll', 'dbsnp', 'Chr', 'Start', 
              'End', 'Tag', 'Author', 'HGMD_Magzine', 'Allname', 'Vol', 
              'Page', 'HGMD_PublishYear', 'HGMD_PMID', 'Reftag', 'Comments', 'ACC_NUM', 
              'New_date',  'Base']
hgmd_use_header = ['Chr', 'Start', 'End', 'HGMD_Disease', 'HGMD_Gene', \
                   'HGMD_Genename', 'HGMD_hgvsAll', 'HGMD_Codon', 'HGMD_PMID',\
                   'HGMD_Magzine', 'HGMD_PublishYear']
#####################################################################
#  table_annovar: protocols, operations
#####################################################################                     
protocols=[refgene, knowngene, ensgene, cytoband, superdup, gwas, snp,\
            cg, ljb, clinvar, cosmic, nci]
operations=['g', 'g', 'g', 'r', 'r', 'r','f', 'f', 'f', 'f', 'f', 'f']
nastring='.'
tblanno_suffix = '_multianno.csv'
tblanno_header=['Chr', 'Start', 'End', 'Ref', 'Alt', \
                'VariantFunction', 'Gene', 'GeneDetail', \
                'ExonicFunction', 'AAChange', 'VariantFunction.k',\
                'Gene.k', 'GeneDetail.k', 'ExonicFunction.k', \
                'AAChange.k', 'VariantFunction.e', 'Gene.e', \
                'GeneDetail.e',    'ExonicFunction.e', 'AAChange.e',\
                'CytoBand', 'GenomicSuperDups', 'GwasCatalog', \
                'dbsnp138',  'CG69', \
                'SIFT_score', 'SIFT_score_converted', 'SIFT_pred', \
                'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred', \
                'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred', \
                'LRT_score', 'LRT_score_converted', 'LRT_pred',  \
                'MutationTaster_score', 'MutationTaster_score_converted', \
                'MutationTaster_pred', \
                'MutationAssessor_score', 'MutationAssessor_score_converted', \
                'MutationAssessor_pred', \
                'FATHMM_score', 'FATHMM_score_converted', 'FATHMM_pred', \
                'RadialSVM_score', 'RadialSVM_score_converted', \
                'RadialSVM_pred', \
                'LR_score', 'LR_pred', 'GERP++', 'PhyloP', 'SiPhy', \
                'clinvar_20140929', 'cosmic70', 'nci60']
#####################################################################
#  merge gene entry info
#####################################################################  
gene_header = ['VariantFunction', 'Gene', 'GeneDetail', 'ExonicFunction', \
               'AAChange', 'VariantFunction.k', 'Gene.k', 'GeneDetail.k', \
               'ExonicFunction.k', 'AAChange.k', 'VariantFunction.e', \
               'Gene.e', 'GeneDetail.e', 'ExonicFunction.e', 'AAChange.e']
region_header = ['CytoBand', 'GenomicSuperDups', 'GwasCatalog']
ljb_pred_header=['SIFT_pred', 'Polyphen2_HDIV_pred', \
                 'Polyphen2_HVAR_pred', 'LRT_pred', \
                 'MutationTaster_pred', 'MutationAssessor_pred', \
                 'FATHMM_pred', 'RadialSVM_pred', 'LR_pred',\
                 'GERP++', 'PhyloP', 'SiPhy']
allele_frequency_header =  g1k_use_header + esp_use_header
other_filter_header = ['clinvar_20140929', 'cosmic70', 'nci60']    
merged_custom_header = basic_header + refgene_header + region_header \
                + ['dbsnp138'] + allele_frequency_header + ljb_pred_header \
                + cadd_header[1:3] + other_filter_header + hgmd_use_header[3:] 
#####################################################################
#  inheritance model filtering
#####################################################################  
denovo_maf=0.001
recessive_maf=0.005
dominant_maf = 0.005
pathogenic_header = ['SIFT_pred', 'Polyphen2_HDIV_pred', \
                        'Polyphen2_HVAR_pred', 'GERP++', 'CADD_phred']
im_header = basic_header + refgene_header + region_header + pathogenic_header \
            + vcf_basic_header 
# other means families missing the proband 
family_type=['singleton', 'duos', 'trios', 'singletonfams', 'duosfams', \
             'triosfams']
allele_num = 6
