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
    input_dataframe = pd.read_csv(input_file, sep="\t", columns = ["CHROM", "POS", "ALT", "REF", "IMP_gt", "TRUTH_gt", "type"])
    frequence = input_dataframe.type.value_counts()
    pct_concordant = frequence['c']/(frequence['c'] + frequence['d'] + frequence['o'])
    pct_missing = frequence['m']/(frequence['c'] + frequence['d'] + frequence['m'])
    pct_only_imputed = frequence['o']/(frequence['c'] + frequence['d'] + frequence['o'])

    result_dict = {'SID':[], 'NC':[frequence['c']], 'ND': [frequence['d']], 'NO': [frequence['o']], 'NM': [frequence['m']], 'PC':[pct_concordant], 'PM':[pct_missing], 'PO':[pct_only_imputed]}
    result_dataframe = pd.DataFrame(result_dict)
    result_dataframe.to_csv(path_out, sep="\t", index = False)

if __name__ == "__main__":
    args = argparser.parse_args()
    input_file = args.in_file_name  
    path_out = args.out_file_name
    sample_ID = args.sample_ID
    calculates_stats(sample_ID, input_file, path_out)
