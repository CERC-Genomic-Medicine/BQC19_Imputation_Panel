#!/usr/bin/env python

import pysam
import argparse
import pandas as pd 
import numpy as np 

argparser = argparse.ArgumentParser(description = 'This script generate analysis based on aggregated all the summary data generated for each evaluated individual.')

argparser.add_argument('-i', '--input_summary_file', metavar = 'file', dest = 'in_file_name', type = str, required = True, help = 'input concatenated summary files.')
argparser.add_argument('-o', '--output_file', metavar = 'file', dest = 'out_file_name', type = str, required = True, help = 'Output file name. Output file will be compressed using gzip.')

def aggregate_all_samples(input_file, path_out):
    all_df = pd.read_csv(input_file, sep="\t", columns=['Sample ID', 'NC', 'ND', 'NO', 'NM', 'PC', 'PM', 'PO'])
    average_concordance = np.mean(all_df['PC'])
    sd_concordance = np.std(all_df['PC'])
    result_dict = {'MEAN' : [average_concordance], 'SD' : [sd_concordance]}
    result_dict.to_csv(path_out, sep = "\t", index = None)

if __name__ == "__main__":
    args = argparser.parse_args()
    input_file = args.in_file_name   
    path_out = args.out_file_name
    aggregate_all_samples(input_file, path_out)
