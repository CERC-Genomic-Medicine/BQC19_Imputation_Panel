#!/usr/bin/env python

import pandas as pd 
import argparse

argparser = argparse.ArgumentParser(description = 
 'This script compares imputation quality analysis of two different reference panel.'
 )

argparser.add_argument('-fq', '--first_quality_file', metavar = 'file', dest = 'in_first_quality_file', type = str, required = True, help = 'Imputation quality file for the first reference panel.')
argparser.add_argument('-sq', '--second_quality_file', metavar = 'file', dest = 'in_second_quality_file', type = str, required = True, help = 'Imputation quality file for the second reference panel.')
argparser.add_argument('-s', '--sample_name', metavar = 'name', dest = 'in_sample_name', type = str, required = True, help = 'name of the sample that analysis performed on.')
argparser.add_argument('-fr', '--first_reference_name', metavar = 'name', dest = 'in_first_reference_name', type = str, required = True, help = 'name of the first reference panel is used for imputation.')
argparser.add_argument('-sr', '--second_reference_name', metavar = 'name', dest = 'in_second_reference_name', type = str, required = True, help = 'name of the second reference panel is used for imputation.')
argparser.add_argument('-c', '--chromosome_name', metavar = 'name', dest = 'in_chr', type = str, required = True, help = 'name of the chromosome is used for processing.')
argparser.add_argument('-m', '--mode', metavar = 'name', dest = 'mode', type = str, required = True, help = 'mode of analysis:{r : raw files, wc : without centromeres, hc : high coverage regions}')

def count_badly_imputed(path):
    df = pd.read_csv(path, sep="\t")
    df_bad = df[df["imputation quality"] == False]
    return len(df_bad)

def count_well_imputed(path):
    df = pd.read_csv(path, sep="\t")
    df_well = df[df["imputation quality"] == True]
    return len(df_well)

def count_shared_variants(path_first, path_second):
    df_first = pd.read_csv(path_first, sep="\t")
    df_second = pd.read_csv(path_second, sep="\t")
    df_merge = df_first.merge(df_second, on=["CHROM", "POS", "REF", "ALT"])
    return len(df_merge)

def count_badly_imputed_first_present_second(path_first, path_second):
    df_first = pd.read_csv(path_first, sep="\t")
    df_first = df_first[df_first["imputation quality"] == False]
    df_second = pd.read_csv(path_second, sep="\t")
    df_merge = df_first.merge(df_second, on=["CHROM", "POS", "REF", "ALT"])
    return len(df_merge)

def count_badly_imputed_in_one_well_imputed_other(path_first, path_second):
    df_first = pd.read_csv(path_first, sep="\t")
    df_second = pd.read_csv(path_second, sep="\t")
    df_first.rename(columns={"imputation quality":"first imputation quality"})
    df_second.rename(columns={"imputation quality":"second imputation quality"})
    df_merge = df_first.merge(df_second)
    return len(df_merge)

if __name__ == "__main__":
    args = argparser.parse_args()
    sample_name = args.in_sample_name   
    path_first = args.in_first_quality_file
    path_second = args.in_second_quality_file
    chrom_name = args.in_chr
    first_reference_name = args.first_reference_name
    second_reference_name = args.second_reference_name
    mode = args.mode

    cbf = count_badly_imputed(path_first)
    cbs = count_badly_imputed(path_second)
    csf = count_well_imputed(path_first)
    css = count_well_imputed(path_second)
    cshv = count_shared_variants(path_first, path_second)
    cbfps = count_badly_imputed_first_present_second(path_first, path_second)
    cbspf = count_badly_imputed_first_present_second(path_second, path_first)
    cwfbs = count_badly_imputed_in_one_well_imputed_other(path_first, path_second)
    cbfws = count_badly_imputed_in_one_well_imputed_other(path_second, path_first)


    result = {"Sample name" : [sample_name], "reference name" : [ref_name],  "chromosome" : [chrom_name], "concordance" : [concordance]}
    df_res = pd.DataFrame(result) 
    df_res.to_csv(sample_name + "_" + chrom_name + "_" + ref_name + "_concordance.txt", sep = "\t", index = None)
    merge_df.to_csv(sample_name + "_" + chrom_name + "_" + ref_name + "_imputation_qualities.txt", sep = "\t", index = None)
