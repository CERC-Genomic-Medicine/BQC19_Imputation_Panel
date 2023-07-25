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

def count_shared_variants(df_first, df_second):
    df_first = df_first[(df_first["type"] == 'WGS_AND_REF_EQ') | (df_first["type"] == 'WGS_AND_REF_GT') | (df_first["type"] == 'WGS_AND_REF_LT')]
    df_second = df_second[(df_second["type"] == 'WGS_AND_REF_EQ') | (df_second["type"] == 'WGS_AND_REF_GT') | (df_second["type"] == 'WGS_AND_REF_LT')]
    df_merge = df_first.merge(df_second, on=["CHROM", "POS", "REF", "ALT"])
    return len(df_merge)

def count_present_in_one_missing_in_the_other(df_first, df_second, first_reference_name, second_reference_name):
    df_first = df_first.rename(columns={"type":"first type"})
    df_second = df_second.rename(columns={"type":"second type"})
    df_merge = df_first.merge(df_second, on=["CHROM", "POS", "REF", "ALT"])
    df_merge = df_merge[((df_merge["first type"] == 'WGS_AND_REF_GT') | (df_merge["first type"] == 'WGS_AND_REF_LT') | (df_merge["first type"] == 'WGS_AND_REF_EQ')) & (df_merge["second type"] == 'WGS')]
    df_merge.to_csv(sample_name + "_present_" + first_reference_name + "_missing_" + second_reference_name + "_pm.txt", sep = "\t", index = None)
    return len(df_merge)

def count_badly_imputed_in_one_well_imputed_other(df_first, df_second, first_reference_name, second_reference_name):
    df_first = df_first.rename(columns={"type":"first type"})
    df_second = df_second.rename(columns={"type":"second type"})
    df_merge = df_first.merge(df_second, on=['CHROM', 'POS', 'ALT', 'REF'])
    df_merge = df_merge[((df_merge["first type"] == 'WGS_AND_REF_GT') | (df_merge["first type"] == 'WGS_AND_REF_LT')) & (df_merge["second type"] == 'WGS_AND_REF_EQ')]
    df_merge.to_csv(sample_name + "_badly_" + first_reference_name + "_well_" + second_reference_name + "_bw.txt", sep = "\t", index = None)
    return len(df_merge)

def count_bad_imputed_in_one_missing_in_other(df_first, df_second, first_reference_name, second_reference_name):
    df_first = df_first.rename(columns={"type":"first type"})
    df_second = df_second.rename(columns={"type":"second type"})
    df_merge = df_first.merge(df_second, on=['CHROM', 'POS', 'ALT', 'REF'])
    df_merge = df_merge[(df_merge["second type"] == 'WGS') & ((df_merge["first type"] == 'WGS_AND_REF_GT') | (df_merge["first type"] == 'WGS_AND_REF_LT'))]
    df_merge.to_csv(sample_name + "_bad_" + first_reference_name + "_missing_" + second_reference_name +  "_bm.txt", sep = "\t", index = None)
    return len(df_merge)

def count_well_imputed_in_one_missing_in_other(df_first, df_second, first_reference_name, second_reference_name):
    df_first = df_first.rename(columns={"type":"first type"})
    df_second = df_second.rename(columns={"type":"second type"})
    df_merge = df_first.merge(df_second, on=['CHROM', 'POS', 'ALT', 'REF'])
    df_merge = df_merge[(df_merge["second type"] == 'WGS') & (df_merge["first type"] == 'WGS_AND_REF_EQ')]
    df_merge.to_csv(sample_name + "_well_" + first_reference_name + "_missing_" + second_reference_name +  "_wm.txt", sep = "\t", index = None)
    return len(df_merge)

def count_well_imputed_both(df_first, df_second):
    df_first = df_first.rename(columns={"type":"first type"})
    df_second = df_second.rename(columns={"type":"second type"})
    df_merge = df_first.merge(df_second, on=['CHROM', 'POS', 'ALT', 'REF'])
    df_merge = df_merge[(df_merge["first type"] == 'WGS_AND_REF_EQ') & (df_merge["second type"] == 'WGS_AND_REF_EQ')]
    df_merge.to_csv(sample_name + "_well_" + first_reference_name + "_well_" + second_reference_name +  "_ww.txt", sep = "\t", index = None)
    return len(df_merge)

def count_badly_imputed_both(df_first, df_second, first_reference_name, second_reference_name):
    df_first = df_first.rename(columns={"type":"first type"})
    df_second = df_second.rename(columns={"type":"second type"})
    df_merge = df_first.merge(df_second, on=['CHROM', 'POS', 'ALT', 'REF'])
    df_merge = df_merge[((df_merge["first type"] == 'WGS_AND_REF_GT') | (df_merge["first type"] == 'WGS_AND_REF_LT')) & ((df_merge["second type"] == 'WGS_AND_REF_GT') | (df_merge["second type"] == 'WGS_AND_REF_LT'))]
    df_merge.to_csv(sample_name + "_badly_" + first_reference_name + "_badly_" + second_reference_name  + "_bb.txt", sep = "\t", index = None)
    return len(df_merge)

def count_missing_both(df_first, df_second, first_reference_name, second_reference_name):
    df_first = df_first.rename(columns={"type":"first type"})
    df_second = df_second.rename(columns={"type":"second type"})
    df_merge = df_first.merge(df_second, on=['CHROM', 'POS', 'ALT', 'REF'])
    df_merge = df_merge[(df_merge["first type"] == 'WGS') & (df_merge["second type"] == 'WGS')]
    df_merge.to_csv(sample_name + "_missing_" + first_reference_name + "_missing_" + second_reference_name  + "_mm.txt", sep = "\t", index = None)
    return len(df_merge)

if __name__ == "__main__":
    args = argparser.parse_args()
    sample_name = args.in_sample_name   
    path_first = args.in_first_quality_file
    path_second = args.in_second_quality_file
    first_reference_name = args.in_first_reference_name
    second_reference_name = args.in_second_reference_name
    

    df_first = pd.read_csv(path_first, sep="\t", usecols=['CHROM', 'POS', 'REF', 'ALT','type'])
    df_second = pd.read_csv(path_second, sep="\t", usecols=['CHROM', 'POS', 'REF', 'ALT', 'type'])

    cshv = count_shared_variants(df_first, df_second)
    cbfws = count_badly_imputed_in_one_well_imputed_other(df_first, df_second, first_reference_name, second_reference_name)
    cwfbs = count_badly_imputed_in_one_well_imputed_other(df_second, df_first, second_reference_name, first_reference_name)
    cwfms = count_well_imputed_in_one_missing_in_other(df_first, df_second, first_reference_name, second_reference_name)
    cwsmf = count_well_imputed_in_one_missing_in_other(df_second, df_first, second_reference_name, first_reference_name)
    cbfms = count_bad_imputed_in_one_missing_in_other(df_first, df_second, first_reference_name, second_reference_name)
    cbsmf = count_bad_imputed_in_one_missing_in_other(df_second, df_first, second_reference_name, first_reference_name)
    cpfms = count_present_in_one_missing_in_the_other(df_first, df_second, first_reference_name, second_reference_name)
    cpsmf = count_present_in_one_missing_in_the_other(df_second, df_first, second_reference_name, first_reference_name)
    cwfws = count_well_imputed_both(df_first, df_second)
    cbfbs = count_badly_imputed_both(df_first, df_second, first_reference_name, second_reference_name)
    cmfms = count_missing_both(df_first, df_second, first_reference_name, second_reference_name)


    result = {"Sample name" : [sample_name], "first reference name" : [first_reference_name], "second reference name" : [second_reference_name],\
      "Shared variants " : [cshv],\
      "Disconcordant " + first_reference_name + " Concordant " + second_reference_name : [cbfws], \
      "Disconcordant " + second_reference_name + " Concordant " + first_reference_name: [cwfbs], \
      "Concordant " + first_reference_name + " Missing " + second_reference_name : [cwfms],\
      "Present " + first_reference_name + " Missing " + second_reference_name : [cpfms],\
      "Present " + second_reference_name + " Missing " + first_reference_name : [cpsmf],\
      "Discordant " + first_reference_name + " Missing " + second_reference_name : [cbfms],\
      "Discordant " + second_reference_name + " Missing " + first_reference_name : [cbsmf],\
      "Concordant " + second_reference_name + " Missing " + first_reference_name : [cwsmf], "Concordant in both" : [cwfws], "Disconcordant in both" : [cbfbs], "Missing in both" : [cmfms]}


    df_res = pd.DataFrame(result) 
    df_res.to_csv(sample_name + "_" + first_reference_name + "_" + second_reference_name +  "_post_imputation_analysis.txt", sep = "\t", index = None)

