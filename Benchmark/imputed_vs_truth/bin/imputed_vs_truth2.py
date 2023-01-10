#!/usr/bin/env python

import pysam
import gzip
import argparse

argparser = argparse.ArgumentParser(description = 'This script compares imputed genotypes against true genotypes.')

argparser.add_argument('-iv', '--imputed_vcf', metavar = 'file', dest = 'in_imp_vcf', type = str, required = True, help = 'VCF file containing imputed data.')
argparser.add_argument('-tv', '--truth_vcf', metavar = 'file', dest = 'in_truth_vcf', type = str, required = True, help = 'VCF file containing truth data.')
argparser.add_argument('-s', '--sample_name', metavar = 'name', dest = 'in_sample_name', type = str, required = True, help = 'Name of the sample that analysis performed on.')
argparser.add_argument('-o', '--output_file', metavar = 'file', dest = 'out_file_name', type = str, required = True, help = 'Output file name. Output file will be compressed using gzip.')



def read_variant(filename, sample_name, imputed_flag, chrom, start, stop):
    with pysam.VariantFile(filename) as ivcf:
        if any(c.startswith('chr') for c in ivcf.header.contigs): # add or remove 'chr' prefix if needed, depending on VCF header
            if not chrom.startswith('chr'):
                chrom = f'chr{chrom}'
        else:
            if chrom.startswith('chr'):
                chrom = chrom[3:]
        ivcf.subset_samples([sample_name])
        for record in ivcf.fetch(contig = chrom, start = start, stop = stop):
            if (imputed_flag == True):
             if not record.info["IMPUTED"]:
                continue
            assert len(record.alts) == 1
            if (len(record.alts) != 1):
                continue
            gt = list(value['GT'] for value in record.samples.values())[0]
            yield (record.pos, record.ref, record.alts[0], gt)


def compare(imputed_gt_filename, truth_gt_filename, sample_name):
    with pysam.VariantFile(imputed_gt_filename) as ivcf:
        chroms = list(ivcf.header.contigs)

    for chrom in chroms:
        imp_variants = read_variant(imputed_gt_filename, sample_name, True, chrom, None, None)
        truth_variants = read_variant(truth_gt_filename, sample_name, False, chrom, None, None)
        imp_variants_buffer = []
        for truth_pos, truth_ref, truth_alt, truth_gt in truth_variants:
            for imp_pos, imp_ref, imp_alt, imp_gt in imp_variants:
                imp_variants_buffer.append((imp_pos, imp_ref, imp_alt, imp_gt))
                if imp_pos > truth_pos:
                    break
            imputed_truth = False 
            while imp_variants_buffer:
                imp_pos, imp_ref, imp_alt, imp_gt = imp_variants_buffer[0]
                if imp_pos < truth_pos:
                    imp_variants_buffer.pop(0)
                    print(chrom, imp_pos, imp_ref, imp_alt, imp_gt) # only imputed
                elif imp_pos == truth_pos:
                    imp_variants_buffer.pop(0)
                    if imp_ref == truth_ref and imp_alt == truth_alt:
                        print(chrom, imp_pos, imp_ref, imp_alt, imp_gt, truth_pos, truth_ref, truth_alt, truth_gt) # imputed and in truth
                        imputed_truth = True
                        break
                    else:
                        print(chrom, imp_pos, imp_ref, imp_alt, imp_gt) # only imputed
                else:
                    break
            if not imputed_truth:
                print( chrom, truth_pos, truth_ref, truth_alt, truth_gt) # only in truth
        
        for imp_pos, imp_ref, imp_alt, imp_gt in imp_variants_buffer:
            print( chrom, imp_pos, imp_ref, imp_alt, imp_gt) # only imputed


if __name__ == "__main__":
    args = argparser.parse_args()
    sample_name = args.in_sample_name   
    path_imputed = args.in_imp_vcf
    path_truth = args.in_truth_vcf

    compare(args.in_imp_vcf, args.in_truth_vcf, sample_name)
 
