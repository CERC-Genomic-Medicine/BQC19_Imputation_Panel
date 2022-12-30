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

def count_badly_imputed(df):
    df_bad = df[df["imputation quality"] == False]
    return len(df_bad)

def count_well_imputed(df):
    df_well = df[df["imputation quality"] == True]
    return len(df_well)

def count_missing(df):
    df_missing = df[df["imputation quality"] == "missing"]
    return len(df_missing)

def count_shared_variants(df_first, df_second):
    df_first = df_first[df_first["imputation quality"] != "missing"]
    df_second = df_second[df_second["imputation quality"] != "missing"]
    df_merge = df_first.merge(df_second, on=["CHROM", "POS", "REF", "ALT"])
    return len(df_merge)

def count_badly_imputed_first_present_second(df_first, df_second):
    df_second = df_second[df_second["imputation quality"] != "missing"]
    df_first = df_first[df_first["imputation quality"] == False]
    df_merge = df_first.merge(df_second, on=["CHROM", "POS", "REF", "ALT"])
    return len(df_merge), df_merge

def count_badly_imputed_in_one_well_imputed_other(df_first, df_second):
    df_first.rename(columns={"imputation quality":"first imputation quality"})
    df_second.rename(columns={"imputation quality":"second imputation quality"})
    df_merge = df_first.merge(df_second)
    df_merge = df_merge[(df_merge["first imputation quality"] == False) and (df_merge["second imputation quality"] == True)]
    return len(df_merge)

def count_present_in_one_missing_in_other(df_first, df_second):
    df_first.rename(columns={"imputation quality":"first imputation quality"})
    df_second.rename(columns={"imputation quality":"second imputation quality"})
    df_merge = df_first.merge(df_second)
    df_merge = df_merge[(df_merge["first imputation quality"] != "missing") and (df_merge["second imputation quality"] == "missing")]
    return len(df_merge)

def count_well_imputed_in_one_missing_in_other(df_first, df_second):
    df_first.rename(columns={"imputation quality":"first imputation quality"})
    df_second.rename(columns={"imputation quality":"second imputation quality"})
    df_merge = df_first.merge(df_second)
    df_merge = df_merge[(df_merge["first imputation quality"] == True) and (df_merge["second imputation quality"] == "missing")]
    return len(df_merge), df_merge

def count_well_imputed_both(df_first, df_second):
    df_first.rename(columns={"imputation quality":"first imputation quality"})
    df_second.rename(columns={"imputation quality":"second imputation quality"})
    df_merge = df_first.merge(df_second)
    df_merge = df_merge[(df_merge["first imputation quality"] == True) and (df_merge["second imputation quality"] == True)]
    return len(df_merge)

def count_badly_imputed_both(df_first, df_second):
    df_first.rename(columns={"imputation quality":"first imputation quality"})
    df_second.rename(columns={"imputation quality":"second imputation quality"})
    df_merge = df_first.merge(df_second)
    df_merge = df_merge[(df_merge["first imputation quality"] == False) and (df_merge["second imputation quality"] == False)]
    return len(df_merge)


def count_missing_both(df_first, df_second):
    df_first.rename(columns={"imputation quality":"first imputation quality"})
    df_second.rename(columns={"imputation quality":"second imputation quality"})
    df_merge = df_first.merge(df_second)
    df_merge = df_merge[(df_merge["first imputation quality"] == "missing") and (df_merge["second imputation quality"] == "missing")]
    return len(df_merge), df_merge

if __name__ == "__main__":
    args = argparser.parse_args()
    sample_name = args.in_sample_name   
    path_first = args.in_first_quality_file
    path_second = args.in_second_quality_file
    chrom_name = args.in_chr
    first_reference_name = args.first_reference_name
    second_reference_name = args.second_reference_name
    mode = args.mode

    df_first = pd.read_csv(path_first, sep="\t")
    df_second = pd.read_csv(path_second, sep="\t")

    cbf = count_badly_imputed(df_first)
    cbs = count_badly_imputed(df_second)
    cwf = count_well_imputed(df_first)
    cws = count_well_imputed(df_second)
    cmf = count_missing(df_first)
    cms = count_missing(df_second)

    cshv = count_shared_variants(df_first, df_second)
    cbfps = count_badly_imputed_first_present_second(df_first, df_second)
    cbspf = count_badly_imputed_first_present_second(df_second, df_first)
    cbfws, df_bfws = count_badly_imputed_in_one_well_imputed_other(df_first, df_second)
    cwfbs, df_wfbs = count_badly_imputed_in_one_well_imputed_other(df_second, df_first)
    cpfms = count_present_in_one_missing_in_other(df_first, df_second)
    cmfps = count_present_in_one_missing_in_other(df_second, df_first)
    cwfms, df_wfms = count_well_imputed_in_one_missing_in_other(df_first, df_second)
    cwsmf, df_wsmf = count_well_imputed_in_one_missing_in_other(df_second, df_first)

    cwfws = count_well_imputed_both(df_first, df_second)
    cbfbs = count_badly_imputed_both(df_first, df_second)
    cmfms, df_mfms = count_missing_both(df_first, df_second)



    result = {"Sample name" : [sample_name], "chromosome" : [chrom_name], "first reference name" : [first_reference_name], "second reference name" : [second_reference_name],\
     "Disconcordant " + first_reference_name : [cbf], "Disconcordant " + second_reference_name : [cbs], "Concordant " + first_reference_name : [cwf], "Concordant " + second_reference_name : [cws],\
     "Missing " + first_reference_name : [cmf], "Missing " + second_reference_name : [cms], "Shared variants " : [cshv], "Disconcordant " + first_reference_name + " Present " + second_reference_name : [cbfps],\
      "Disconcordant " + second_reference_name + " Present " + first_reference_name: [cbspf], "Disconcordant " + first_reference_name + " Concordant " + second_reference_name : [cbfws], \
      "Disconcordant " + second_reference_name + " Concordant " + first_reference_name: [cwfbs], "Present " + first_reference_name + " Missing " + second_reference_name : [cpfms], \
      "Present " + first_reference_name + " Missing " + second_reference_name: [cmfsp], "Concordant " + first_reference_name + " Missing " + second_reference_name : [cwfms],\
      "Concordant " + second_reference_name + " Missing " + first_reference_name : [cwsmf], "Concordant in both" : [cwfws], "Disconcordant in both" : [cbfbs], "Missing in both" : [cmfms]}


    df_res = pd.DataFrame(result) 
    df_res.to_csv(sample_name + "_" + chrom_name + "_" + first_reference_name + "_" second_reference_name + "_post_imputation_analysis.txt", sep = "\t", index = None)
    df_bfws.to_csv(sample_name + "_" + chrom_name + "_badly_" + first_reference_name + "_well_" + second_reference_name + ".txt", sep = "\t", index = None)
    df_wfbs.to_csv(sample_name + "_" + chrom_name + "_well_" + first_reference_name + "_badly_" + second_reference_name + ".txt", sep = "\t", index = None)
    df_wfms.to_csv(sample_name + "_" + chrom_name + "_well_" + first_reference_name + "_missing_" + second_reference_name + ".txt", sep = "\t", index = None)
    df_wsmf.to_csv(sample_name + "_" + chrom_name + "_missing_" + first_reference_name + "_well_" + second_reference_name + ".txt", sep = "\t", index = None)
    df_mfms.to_csv(sample_name + "_" + chrom_name + "_missing_" + first_reference_name + "_missing_" + second_reference_name + ".txt", sep = "\t", index = None)
