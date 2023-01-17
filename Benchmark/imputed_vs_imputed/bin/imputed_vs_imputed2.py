#!/usr/bin/env python

import pysam
import argparse

argparser = argparse.ArgumentParser(description = 
 'This script compares imputation quality analysis of two different reference panel.'
 )

argparser.add_argument('-fq', '--first_quality_file', metavar = 'file', dest = 'in_first_quality_file', type = str, required = True, help = 'Imputation quality file for the first reference panel.')
argparser.add_argument('-sq', '--second_quality_file', metavar = 'file', dest = 'in_second_quality_file', type = str, required = True, help = 'Imputation quality file for the second reference panel.')
argparser.add_argument('-s', '--sample_name', metavar = 'name', dest = 'in_sample_name', type = str, required = True, help = 'name of the sample that analysis performed on.')
argparser.add_argument('-fr', '--first_reference_name', metavar = 'name', dest = 'in_first_reference_name', type = str, required = True, help = 'name of the first reference panel is used for imputation.')
argparser.add_argument('-sr', '--second_reference_name', metavar = 'name', dest = 'in_second_reference_name', type = str, required = True, help = 'name of the second reference panel is used for imputation.')
argparser.add_argument('-m', '--mode', metavar = 'name', dest = 'mode', type = str, required = True, help = 'mode of analysis:{r : raw files, wc : without centromeres, hc : high coverage regions}')



def read_variant(filename):
    with gzip.open(filename, 'rt') as fr:
        yield (fr.readline().rstrip().split())


def compare(first_reference_qt_filename, second_reference_qt_filename, path_out):
        first_variants = read_variant(first_reference_qt_filename)
        second_variants = read_variant(second_reference_qt_filename)
        second_variants_buffer = []
        only_present_in_second_ref = 0
        only_present_in_first_ref = 0
        both_missing = 0
        first_missing_second_concordant = 0
        first_missing_second_discordant = 0
        first_missing_second_missing = 0
        first_concordant_second_missing = 0
        first_discordant_second_missing = 0
        first_concordant_second_concordant = 0
        first_concordant_second_discordant = 0
        first_discordant_second_discordant = 0
        with pysam.BGZFile(path_out, 'w')  as fw:
            for f_chrom, f_pos, f_ref, f_alt, f_imp_gt, f_truth_gt, f_qt in first_variants:
                for s_chrom, s_pos, s_ref, s_alt, s_imp_gt, s_truth_gt, s_qt in second_variants:
                    second_variants_buffer.append((s_chrom, s_pos, s_ref, s_alt, s_imp_gt, s_truth_gt, s_qt))
                    if(f_chrom != s_chrom):
                        raise Exception()
                    if s_pos > f_pos:
                        break
                
                present_both = False
                while second_variants_buffer:
                    s_chrom, s_pos, s_ref, s_alt, s_imp_gt, s_truth_gt, s_qt = second_variants_buffer[0]
                    if s_pos < f_pos:
                        second_variants_buffer.pop(0)
                        only_present_in_second_ref += 1
                    elif s_pos == f_pos:
                        second_variants_buffer.pop(0)
                        if f_ref == s_ref and f_alt == s_alt:
                            present_both = True
                            if(f_qt == "only truth"):
                                if(s_qt == "only truth"):
                                    first_missing_second_missing += 1
                                if(s_qt == "concordant"):
                                    first_missing_second_concordant +=1
                                if(s_qt == "discordant"):
                                    first_missing_second_discordant += 1
                            if(s_qt == "only truth"):
                                if(f_qt == "concordant"):
                                    first_concordant_second_missing += 1
                                if(f_qt == "disconcordant"):
                                    first_discordant_second_missing += 1
                        
                            if(f_qt == "concordant"):
                                if(s_qt == "concordant"):
                                    first_concordant_second_concordant += 1
                                if(s_qt == "diconcordant"):
                                    first_concordant_second_discordant += 1
                            if(s_qt == "disconcordant"):
                                if(f_qt == "disconcordant"):
                                    first_discordant_second_discordant += 1

                        else:
                            only_present_in_second_ref += 1
                    else:
                        break    
                if (present_both == False):
                    only_present_in_first_ref += 1
            for s_chrom, s_pos, s_ref, s_alt, s_imp_gt, s_truth_gt, s_qt in second_variants:
                    second_variants_buffer.append((s_chrom, s_pos, s_ref, s_alt, s_imp_gt, s_truth_gt, s_qt))  
            for s_chrom, s_pos, s_ref, s_alt, s_imp_gt, s_truth_gt, s_qt in second_variants_buffer:
                    only_present_in_second_ref += 1

if __name__ == "__main__":
    #to be completed
 
