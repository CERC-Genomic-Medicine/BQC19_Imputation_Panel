import pysam
import pandas as pd
import sys
import argparse

argparser = argparse.ArgumentParser(description = 
 'This script extracts OmniExpress positions.'
 )
argparser.add_argument('-p', '--positions_array', metavar = 'file', dest = 'in_position_txt', type = str, required = True, help = 'file containing position of OmniExpress Array.')
argparser.add_argument('-i', '--input_vcf_file', metavar = 'file', dest = 'in_vcf_file', type = str, required = True, help = 'input vcf file.')
argparser.add_argument('-o', '--output_vcf_file', metavar = 'name', dest = 'out_vcf_file', type = str, required = True, help = 'output vcf file.')


def read_by_sample(file_name):
    with pysam.VariantFile(file_name) as vcf:
        for record in vcf:
            yield {'CHROM':record.chrom, "POS": record.pos, "REF":record.ref, "ALT": record.alts[0]}

if __name__ == '__main__':
    args = argparser.parse_args()
    input_file = args.in_vcf_file
    file_prep_array = args.out_vcf_file
    df_position = pd.read_csv(args.in_position_txt, sep ='\t', names = ["CHROM", "POS"])
    df_unrelated = pd.DataFrame(read_by_sample(input_file))
    df_merge = df_unrelated.merge(df_position)
    with pysam.VariantFile(input_file, 'r') as unrelated_vcf:
            vcf_out = pysam.VariantFile(file_prep_array, 'w', header=unrelated_vcf.header)
            for record in unrelated_vcf:
                df_rec = df_merge[df_merge["POS"] == record.pos]
                if(not df_rec.empty):
                    vcf_out.write(record)
