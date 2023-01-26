#!/usr/bin/env python

import pysam
import argparse
import pandas as pd 

argparser = argparse.ArgumentParser(description = 'This script calculates some stats for each sample based on the comparison between imputed and truth files.')

argparser.add_argument('-i', '--input_file', metavar = 'file', dest = 'in_file_name', type = str, required = True, help = 'text file containing the labels for each variant, whether it is\
 well imputed, badly imputed, missing or only present in the reference panel but in truth VCFs.')
argparser.add_argument('-s', '--sample_ID', metavar = 'name', dest = 'in_sample_ID', type = str, required = True, help = 'ID of the sample that is used for evaluation.')
argparser.add_argument('-o', '--output_file', metavar = 'file', dest = 'out_file_name', type = str, required = True, help = 'Output file which contains summary of numbers of different types of variant.')

def calculates_stats(sample_ID, input_file, path_out):
    input_dataframe = pd.read_csv(input_file, sep="\t", header = None)
    input_dataframe.columns = ["CHROM", "POS", "ALT", "REF", "IMP_gt", "TRUTH_gt", "tag", "zero_or_not", "type"]
    
    counts_type = input_dataframe.type.value_counts()
    list_index = list(counts_type.index)
    
    Only_REF = counts_type['Only_REF'] if ('Only_REF' in list_index) else  0
    WGS_AND_REF = counts_type['WGS_AND_REF'] if ('WGS_AND_REF' in list_index) else  0
    Only_WGS = counts_type['Only_WGS'] if ('Only_WGS' in list_index) else  0

    counts_tag = input_dataframe.tag.value_counts()
    list_index = list(counts_tag.index)

    WGS_AND_REF_EQ = counts_tag['WGS_AND_REF_EQ'] if ('WGS_AND_REF_EQ' in list_index) else 0
    WGS_AND_REF_LT = counts_tag['WGS_AND_REF_LT'] if ('WGS_AND_REF_LT' in list_index) else 0
    WGS_AND_REF_GT = counts_tag['WGS_AND_REF_GT'] if ('WGS_AND_REF_GT' in list_index) else 0

    counts_zero_or_not = input_dataframe.zero_or_not.value_counts()
    list_index = list(counts_zero_or_not.index)

    WGS_Zero = counts_zero_or_not['ZGT'] if ('ZGT' in list_index) else 0
    WGS_Non_Zero = counts_zero_or_not['NZGT'] if ('NZGT' in list_index) else 0
    IMP_Zero = counts_zero_or_not['ZIMP'] if ('ZIMP' in list_index) else 0
    IMP_Non_Zero = counts_zero_or_not['NZIMP'] if ('NZIMP' in list_index) else 0

    result_dict = {'Sample ID':[sample_ID], 'WGS':[Only_WGS], 'REF': [Only_REF], 'WGS_AND_REF': [WGS_AND_REF],\
     'WGS_AND_REF_EQ': [WGS_AND_REF_EQ], 'WGS_AND_REF_LT':[WGS_AND_REF_LT], 'WGS_AND_REF_GT':[WGS_AND_REF_GT],\
     'REF_0ALT' : [IMP_Zero], 'REF_WALT' : [IMP_Non_Zero], 'WGS_0ALT':[WGS_Zero], 'WGS_Non_Zero':[WGS_Non_Zero]}
    result_dataframe = pd.DataFrame(result_dict)
    result_dataframe.to_csv(path_out, sep="\t", index = False)

if __name__ == "__main__":
    args = argparser.parse_args()
    input_file = args.in_file_name  
    path_out = args.out_file_name
    sample_ID = args.in_sample_ID
    calculates_stats(sample_ID, input_file, path_out)
