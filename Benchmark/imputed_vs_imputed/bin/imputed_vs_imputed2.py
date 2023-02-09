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



def read_variant(filename, chr):
    with gzip.open(filename, 'rt') as fr:
        line = fr.readline().rstrip().split()
        if(line[0] == chr):
            yield (line)

def compare(sample_name, first_reference_qt_filename, second_reference_qt_filename, path_present_first_missing_second, path_well_first_badly_second, path_badly_first_well_second, path_res, first_reference_name, second_reference_name):
    
    first_present_second_missing = 0
    first_concordant_second_discordant = 0
    first_discordant_second_discordant = 0
    for chr in ['chr'+str(i) for i in range(1, 23)]:
        first_variants = read_variant(first_reference_qt_filename, chr)
        second_variants = read_variant(second_reference_qt_filename, chr)
        second_variants_buffer = []
        with pysam.VariantFile(path_present_first_missing_second, 'w') as pfms, pysam.VariantFile(path_well_first_badly_second, 'w') as wfbs, pysam.VariantFile(path_badly_first_well_second, 'w') as bfws:
            pfms.write('##fileformat=VCFv4.1\n')
            pfms.write(f'##ref1={first_reference_name}\n')
            pfms.write(f'##ref2={second_reference_name}\n')
            pfms.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            wfbs.write('##fileformat=VCFv4.1\n')
            wfbs.write(f'##ref1={first_reference_name}\n')
            wfbs.write(f'##ref2={second_reference_name}\n')
            wfbs.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            bfws.write('##fileformat=VCFv4.1\n')
            bfws.write(f'##ref1={first_reference_name}\n')
            bfws.write(f'##ref2={second_reference_name}\n')
            bfws.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            present_in_first_ref = 0
            well_first_bad_second = 0
            bad_first_well_second = 0
            for f_chrom, f_pos, f_ref, f_alt, f_imp_gt, f_truth_gt, f_qt in first_variants:
                for s_chrom, s_pos, s_ref, s_alt, s_imp_gt, s_truth_gt, s_qt in second_variants:
                    second_variants_buffer.append((s_chrom, s_pos, s_ref, s_alt, s_imp_gt, s_truth_gt, s_qt))
                    if(f_chrom != s_chrom):
                        raise Exception()
                    if s_pos > f_pos:
                        break
                
                while second_variants_buffer:
                    s_chrom, s_pos, s_ref, s_alt, s_imp_gt, s_truth_gt, s_qt = second_variants_buffer[0]
                    if s_pos < f_pos:
                        second_variants_buffer.pop(0)

                    elif s_pos == f_pos:
                        second_variants_buffer.pop(0)
                        if f_ref == s_ref and f_alt == s_alt:
                            if((f_qt == "WGS_AND_REF_EQ") or (f_qt == "WGS_AND_REF_GT") or (f_qt == "WGS_AND_REF_LT")):
                                if((s_qt != "WGS_AND_REF_EQ") and (s_qt != "WGS_AND_REF_GT") and (s_qt != "WGS_AND_REF_LT")):
                                    present_in_first_ref +=1
                                    pfms.write(f'{s_chrom}\t{f_pos}\t.\t{f_ref}\t{f_alt}\t.\t.\t.\n')

                            if((f_qt == "WGS_AND_REF_EQ")):
                                if((s_qt == "WGS_AND_REF_GT") or (s_qt == "WGS_AND_REF_LT")):
                                    well_first_bad_second += 1
                                    wfbs.write(f'{s_chrom}\t{f_pos}\t.\t{f_ref}\t{f_alt}\t.\t.\t.\n')

                            if((s_qt == "WGS_AND_REF_EQ")):
                                if((f_qt == "WGS_AND_REF_GT") or (f_qt == "WGS_AND_REF_LT")):
                                    bad_first_well_second += 1
                                    bfws.write(f'{s_chrom}\t{f_pos}\t.\t{f_ref}\t{f_alt}\t.\t.\t.\n')

                    else:
                        break 
        
    
    result = {'sample ID':[sample_id], 'present_'+first_reference_name+'_missing_'+second_reference_name:[present_in_first_ref], 'well_'+first_reference_name+'_bad_'+second_reference_name:[well_first_bad_second], 'well_'+second_reference_name+'_bad_'+first_reference_name:[bad_first_well_second]}
    df_re = pd.DataFrame(result)
    df_re.to_csv(path_res, sep='\t', index=None)

if __name__ == "__main__":
    args = argparser.parse_args()
    sample_name = args.in_sample_name   
    path_first_reference_file = args.in_first_quality_file
    path_second_reference_file = args.in_second_quality_file
    first_reference_name = args.in_first_reference_name
    second_reference_name = args.in_second_reference_name
    path_present_first_missing_second = sample_name + '_present_' + first_reference_name + '_missing_' + second_reference_name + '_pm.vcf.gz'
    path_well_first_badly_second = sample_name + '_well_' + first_reference_name + '_badly_' + second_reference_name + '_wb.vcf.gz'
    path_badly_first_well_second = sample_name + '_badly_' + first_reference_name + '_well_' + second_reference_name + '_bw.vcf.gz'
    path_res = sample_name + '_post_analysis.txt'
    compare(sample_name, path_first_reference_file, path_second_reference_file, path_present_first_missing_second, path_well_first_badly_second, path_badly_first_well_second, first_reference_name, second_reference_name, path_res)
