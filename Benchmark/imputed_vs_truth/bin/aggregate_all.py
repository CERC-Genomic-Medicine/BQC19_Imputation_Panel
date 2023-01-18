#!/usr/bin/env python

import pysam
import argparse

argparser = argparse.ArgumentParser(description = 'This script compares imputed genotypes against true genotypes.')

argparser.add_argument('-iv', '--imputed_vcf', metavar = 'file', dest = 'in_imp_vcf', type = str, required = True, help = 'VCF file containing imputed data.')
argparser.add_argument('-tv', '--truth_vcf', metavar = 'file', dest = 'in_truth_vcf', type = str, required = True, help = 'VCF file containing truth data.')
argparser.add_argument('-s', '--sample_name', metavar = 'name', dest = 'in_sample_name', type = str, required = True, help = 'Name of the sample that analysis performed on.')
argparser.add_argument('-o', '--output_file', metavar = 'file', dest = 'out_file_name', type = str, required = True, help = 'Output file name. Output file will be compressed using gzip.')

def aggregate

if __name__ == "__main__":
    args = argparser.parse_args()
    sample_name = args.in_sample_name   
    path_out = args.out_file_name
