#!/usr/bin/env python

import pysam
import glob
import argparse
import gzip 

argparser = argparse.ArgumentParser(description = 
 'Count the number of each variant per ancestry group.')

argparser.add_argument('-p', '--path_to_all', metavar = 'file', dest = 'in_first_quality_file', type = str, required = True, help = 'path to all vcf files for individuals.')
argparser.add_argument('-s', '--sample_name', metavar = 'name', dest = 'in_sample_name', type = str, required = True, help = 'name of the sample that analysis performed on.')
argparser.add_argument('-fr', '--first_reference_name', metavar = 'name', dest = 'in_first_reference_name', type = str, required = True, help = 'name of the first reference panel is used for imputation.')
argparser.add_argument('-sr', '--second_reference_name', metavar = 'name', dest = 'in_second_reference_name', type = str, required = True, help = 'name of the second reference panel is used for imputation.')

def load_vcf(self):
        with pysam.VariantFile(self.path) as vcf:
            vcf.subset_samples([self.sample_ID])
            for record in vcf:
                assert len(record.alts) == 1
                if(len(record.alts) != 1):
                    continue
                if(len(record.ref) != 1 or len(record.alts[0]) != 1):
                    continue
                chrom = record.chrom if record.chrom.startswith("chr") else "chr" + record.chrom
                gt = list(value['GT'] for value in record.samples.values())[0]
                yield {'CHROM':chrom, "POS": record.pos, "REF":record.ref, "ALT":record.alts[0]}

def count_number_of_each_variant(path_to_all, out_path):
    list_df_samples = []
    for file in glob.glob(path_to_all):
        df_sample = pd.DataFrame(load_vcf)
        list_df_samples.append(df_sample)
    df_concat = pd.concat(df_list).groupby(["CHROM", "POS", "REF", "ALT"]).size().reset_index(name = "count")
    df_concat = df_concat.sort_values(by=['count'],  ascending=False)
    with gzip.open(out_path, 'wt') as union_set:
            union_set.write('##fileformat=VCFv4.1\n')
            union_set.write('##INFO=<ID=VC,Number=1,Type=Float,Description="Count of each variant in subset samples">\n')
            union_set.write(f'##ref1={first_reference_name}\n')
            union_set.write(f'##ref2={second_reference_name}\n')
            union_set.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            for index, row in df_concat.iterrows():
                union_set.write(f'{row['CHROM']}\t{row['POS']}\t.\t{row['REF']}\t{row['ALT']}\t.\t.\tVC={row['count']}\n')

