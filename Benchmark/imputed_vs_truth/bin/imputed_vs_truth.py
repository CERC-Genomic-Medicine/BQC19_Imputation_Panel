#!/usr/bin/env python

import pysam
import pandas as pd 
import os
import argparse

argparser = argparse.ArgumentParser(description = 
 'This script compares imputed data with truth data and calculate imputation concordance based on these analysis'
 )

argparser.add_argument('-iv', '--imputed_vcf', metavar = 'file', dest = 'in_imp_vcf', type = str, required = True, help = 'VCF file containing imputed data.')
argparser.add_argument('-tv', '--truth_vcf', metavar = 'file', dest = 'in_truth_vcf', type = str, required = True, help = 'VCF file containing truth data.')
argparser.add_argument('-s', '--sample_name', metavar = 'name', dest = 'in_sample_name', type = str, required = True, help = 'name of the sample that analysis performed on.')
argparser.add_argument('-r', '--reference_name', metavar = 'name', dest = 'in_reference_name', type = str, required = True, help = 'name of the reference panel is used for imputation.')
argparser.add_argument('-c', '--chromosome_name', metavar = 'name', dest = 'in_chr', type = str, required = True, help = 'name of the chromosome is used for processing.')


def load_imputed(path, sample_ID, imputed_flag = True):
        with pysam.VariantFile(path) as vcf:
            vcf.subset_samples([sample_ID])
            for record in vcf:
                if(imputed_flag == True):
                    if (not record.info["IMPUTED"]):
                        continue
                assert len(record.alts) == 1
                if(len(record.alts) != 1):
                    continue
                if(len(record.ref) != 1 or len(record.alts[0]) != 1):
                    continue
                chrom = record.chrom if record.chrom.startswith("chr") else "chr" + record.chrom
                gt = list(value['GT'] for value in record.samples.values())[0]
                yield {'CHROM':chrom, "POS": record.pos, "REF":record.ref, "ALT":\
                           record.alts[0], "GT_Imputed":gt}


def load_truth(path, sample_ID):
        with pysam.VariantFile(path) as vcf:
            vcf.subset_samples([sample_ID])
            for record in vcf:
                assert len(record.alts) == 1
                if(len(record.alts) != 1):
                    continue
                if(len(record.ref) != 1 or len(record.alts[0]) != 1):
                    continue
                gt = list(value['GT'] for value in record.samples.values())[0]
                yield {'CHROM':record.chrom, "POS": record.pos, "REF":record.ref, "ALT": \
                       record.alts[0], "GT_truth":gt}


def ImputedVsTruth(path_truth, path_imputed, sample_ID, imputed_flag = True):
    imputed_df = pd.DataFrame(load_imputed(path_imputed, sample_ID, imputed_flag))
    truth_df = pd.DataFrame(load_truth(path_truth, sample_ID))
    merge_df = imputed_df.merge(truth_df, on = ["CHROM", "POS", "REF", "ALT"])
    merge_df["# alt allele truth"] = pd.Series(x.count(1) for x in merge_df["GT_truth"])
    merge_df["# alt allele imputed"] = pd.Series(x.count(1) for x in merge_df["GT_Imputed"])
    correct_count = list(merge_df["# alt allele truth"] == merge_df["# alt allele imputed"]).count(True)
    concordance = (correct_count/len(merge_df))*100
    merge_df["imputation quality"] = pd.Series(merge_df["# alt allele truth"] == merge_df["# alt allele imputed"])
    merge_df = merge_df[["CHROM", "POS", "REF", "ALT", "imputation quality"]]
    truth_only_df = imputed_df.merge(truth_df, on = ["CHROM", "POS", "REF", "ALT"], indicator = True)
    truth_only_df = truth_only_df[truth_only_df["_merge"] == "right_only"]
    truth_only_df = truth_only_df[["CHROM", "POS", "REF", "ALT"]]
    truth_only_df["imputation quality"] = pd.Series(['missing' for i in range(len(truth_only_df))])
    final_df = pd.concat([merge_df, truth_only_df])
    final_df = final_df.sort_values(by=["POS"])
    missing = (len(truth_only_df)/len(truth_df))*100
    return missing, concordance, final_df



if __name__ == "__main__":
    args = argparser.parse_args()
    sample_name = args.in_sample_name   
    path_imputed = args.in_imp_vcf
    path_truth = args.in_truth_vcf
    chrom_name = args.in_chr
    ref_name = args.in_reference_name
    missing, concordance, merge_df = ImputedVsTruth(path_truth, path_imputed, sample_name)
    result = {"Sample name" : [sample_name], "reference name" : [ref_name],  "chromosome" : [chrom_name], "concordance" : [concordance], "missing" : [missing]}
    df_res = pd.DataFrame(result) 
    df_res.to_csv(sample_name + "_" + chrom_name + "_" + ref_name + "_concordance.txt", sep = "\t", index = None)
    merge_df.to_csv(sample_name + "_" + chrom_name + "_" + ref_name + "_imputation_qualities.txt", sep = "\t", index = None)
