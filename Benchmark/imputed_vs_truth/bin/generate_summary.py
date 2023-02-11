#!/usr/bin/env python

import pysam
import argparse
import pandas as pd 

argparser = argparse.ArgumentParser(description = 'This script calculates some stats for each sample based on the comparison between imputed and truth files.')

argparser.add_argument('-i', '--input_file', metavar = 'file', dest = 'in_file_name', type = str, required = True, help = 'text file containing the labels for each variant. labels are in [REF, WGS, WGS_0ALT, \
REF_0ALT, WGS_0ALT_AND_REF_EQ, WGS_0ALT_AND_REF_LT, WGS_AND_REF_EQ, WGS_AND_REF_LT]')
argparser.add_argument('-s', '--sample_ID', metavar = 'name', dest = 'in_sample_ID', type = str, required = True, help = 'ID of the sample that is used for evaluation.')
argparser.add_argument('-o', '--output_file', metavar = 'file', dest = 'out_file_name', type = str, required = True, help = 'Output file which contains summary of numbers of different types of variant and corcondance results.')

def calculates_stats(sample_ID, input_file, path_out):
    input_dataframe = pd.read_csv(input_file, sep="\t", usecols = ['type'])
    
    counts_type = input_dataframe.type.value_counts()
    list_index = list(counts_type.index)
    
    Only_REF = counts_type['REF'] if ('REF' in list_index) else  0
    Only_WGS = counts_type['WGS'] if ('WGS' in list_index) else  0
    Only_WGS_0ALT = counts_type['WGS_0ALT'] if ('WGS_0ALT' in list_index) else  0
    Only_REF_0ALT = counts_type['REF_0ALT'] if ('REF_0ALT' in list_index) else  0
    WGS_0ALT_AND_REF_EQ = counts_type['WGS_0ALT_AND_REF_EQ'] if ('WGS_0ALT_AND_REF_EQ' in list_index) else  0
    WGS_0ALT_AND_REF_LT = counts_type['WGS_0ALT_AND_REF_LT'] if ('WGS_0ALT_AND_REF_LT' in list_index) else  0
    WGS_AND_REF_EQ = counts_type['WGS_AND_REF_EQ'] if ('WGS_AND_REF_EQ' in list_index) else 0
    WGS_AND_REF_LT = counts_type['WGS_AND_REF_LT'] if ('WGS_AND_REF_LT' in list_index) else 0
    WGS_AND_REF_GT = counts_type['WGS_AND_REF_GT'] if ('WGS_AND_REF_GT' in list_index) else 0

    WGS_AND_REF = WGS_AND_REF_EQ + WGS_AND_REF_LT + WGS_AND_REF_GT
    WGS = (Only_WGS + WGS_AND_REF_EQ + WGS_AND_REF_LT + WGS_AND_REF_GT)
    AA_Concordance_FRAC = (WGS_AND_REF_EQ / WGS_AND_REF)
    AA_Concordance = WGS_AND_REF_EQ
    Coverage = abs(WGS - WGS_AND_REF)/WGS
    REF = (Only_REF + WGS_0ALT_AND_REF_EQ + WGS_0ALT_AND_REF_LT + Only_REF_0ALT)
    RA_Discordance_FRAC = (Only_REF + WGS_0ALT_AND_REF_LT) / REF
    RA_Discordance = Only_REF + WGS_0ALT_AND_REF_LT

    result_dict = {'Sample ID':[sample_ID], 'WGS':[WGS_AND_REF_EQ + WGS_AND_REF_GT + WGS_AND_REF_LT + WGS], 'REF': [REF],\
    'WGS_AND_REF' : [WGS_AND_REF], 'WGS_AND_REF_EQ': [WGS_AND_REF_EQ],\
    'WGS_AND_REF_LT':[WGS_AND_REF_LT], 'WGS_AND_REF_GT':[WGS_AND_REF_GT],'REF_0ALT' : [Only_REF_0ALT],\
    'WGS_0ALT':[Only_WGS_0ALT + WGS_0ALT_AND_REF_EQ + WGS_0ALT_AND_REF_LT], 'AA_Concordance_FRAC':[AA_Concordance_FRAC],\
    'AA_Concordance':[AA_Concordance], 'Coverage':[Coverage], 'RA_Discordance_FRAC':[RA_Discordance_FRAC], 'RA_Discordance':[RA_Discordance]}

    result_dataframe = pd.DataFrame(result_dict)
    result_dataframe.to_csv(path_out, sep="\t", index = False)

if __name__ == "__main__":
    args = argparser.parse_args()
    input_file = args.in_file_name  
    path_out = args.out_file_name
    sample_ID = args.in_sample_ID
    calculates_stats(sample_ID, input_file, path_out)
