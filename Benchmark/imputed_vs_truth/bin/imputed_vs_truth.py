#!/usr/bin/env python

import pysam
import argparse

argparser = argparse.ArgumentParser(description = 'This script compares imputed genotypes against true genotypes.')

argparser.add_argument('-iv', '--imputed_vcf', metavar = 'file', dest = 'in_imp_vcf', type = str, required = True, help = 'VCF file containing imputed data.')
argparser.add_argument('-tv', '--truth_vcf', metavar = 'file', dest = 'in_truth_vcf', type = str, required = True, help = 'VCF file containing truth data.')
argparser.add_argument('-s', '--sample_name', metavar = 'name', dest = 'in_sample_name', type = str, required = True, help = 'Name of the sample that analysis performed on.')
argparser.add_argument('-o', '--output_file', metavar = 'file', dest = 'out_file_name', type = str, required = True, help = 'Output file name. Output file will be compressed using gzip.')


def read_variant(filename, sample_name, chrom, start, stop):
    with pysam.VariantFile(filename) as ivcf:
        if any(c.startswith('chr') for c in ivcf.header.contigs): # add or remove 'chr' prefix if needed, depending on VCF header
            if not chrom.startswith('chr'):
                chrom = f'chr{chrom}'
        else:
            if chrom.startswith('chr'):
                chrom = chrom[3:]
        ivcf.subset_samples([sample_name])
        for record in ivcf.fetch(contig = chrom, start = start, stop = stop):
            assert len(record.alts) == 1
            if (len(record.alts) != 1):
                continue
            gt = (list(value['GT'] for value in record.samples.values())[0]).count(1)
            yield (record.pos, record.ref, record.alts[0], gt)



def compare(imputed_gt_filename, truth_gt_filename, sample_name, path_out):
    with pysam.VariantFile(imputed_gt_filename) as ivcf:
        chroms = list(ivcf.header.contigs)

    with pysam.BGZFile(path_out, 'w')  as fw:
        fw.write((f"CHROM\tPOS\tREF\tALT\tIMP_gt\tTRUTH_gt\ttype\n").encode())
        for chrom in chroms:
                imp_variants = read_variant(imputed_gt_filename, sample_name, chrom, None, None)
                truth_variants = read_variant(truth_gt_filename, sample_name, chrom, None, None)
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
                            if(imp_gt != 0):
                                fw.write((f"{chrom}\t{imp_pos}\t{imp_ref}\t{imp_alt}\t{imp_gt}\t{None}\tREF\n").encode()) # variant is only present in imputed files, non-zero genotype [(1, 0), (1, 1), (0, 1)].
                            else:
                                fw.write((f"{chrom}\t{imp_pos}\t{imp_ref}\t{imp_alt}\t{imp_gt}\t{None}\tREF_0ALT\n").encode()) # variant is only present in imputed files, zero genotype (0, 0).
                        elif imp_pos == truth_pos:
                            imp_variants_buffer.pop(0)
                            if ((imp_ref == truth_ref) and (imp_alt == truth_alt)):
                                imputed_truth = True
                                if(imp_gt == truth_gt):
                                        if(truth_gt != 0):
                                            fw.write((f"{chrom}\t{imp_pos}\t{imp_ref}\t{imp_alt}\t{imp_gt}\t{truth_gt}\tWGS_AND_REF_EQ\n").encode()) # vartiant is present in both imputed and truth files, non-zero genotype [(1, 0), (1, 1), (0, 1)], number of ALT alleles is equal in imputed vs truth.
                                            break
                                        else:
                                            fw.write((f"{chrom}\t{imp_pos}\t{imp_ref}\t{imp_alt}\t{imp_gt}\t{truth_gt}\tWGS_0ALT_AND_REF_EQ\n").encode()) # variant is present in both imputed and truth files, zero genotype (0, 0), number of ALT alleles is equal in imputed vs truth which is equal to zero in this case.
                                            break
                                elif(truth_gt < imp_gt):
                                        if(truth_gt != 0):
                                            fw.write((f"{chrom}\t{imp_pos}\t{imp_ref}\t{imp_alt}\t{imp_gt}\t{truth_gt}\tWGS_AND_REF_LT\n").encode()) # variant is present in both imputed and truth files, non-zero genotype [(1, 0), (1, 1), (0, 1)] for truth, number of ALT alleles in truth is less than number of ALT alleles in imputed files.
                                            break
                                        else:
                                            fw.write((f"{chrom}\t{imp_pos}\t{imp_ref}\t{imp_alt}\t{imp_gt}\t{truth_gt}\tWGS_0ALT_AND_REF_LT\n").encode()) # variant is present in both imputed and truth files, zero genotype [(0, 0)] for truth, number of ALT alleles in truth is less than number of ALT alleles in imputed files (this means that number of alth allels in imputed files for this variant is non-zero)
                                            break
                                else:
                                        fw.write((f"{chrom}\t{imp_pos}\t{imp_ref}\t{imp_alt}\t{imp_gt}\t{truth_gt}\tWGS_AND_REF_GT\n").encode()) # variant is present in both imputed and truth files, non-zero genotype [(1, 0), (1, 1), (0, 1)] for truth, number of ALT alleles in truth is greater than number of ALT alleles in imputed files.
                                        break
                            
                            else:
                                if(imp_gt != 0):
                                    fw.write((f"{chrom}\t{imp_pos}\t{imp_ref}\t{imp_alt}\t{imp_gt}\t{None}\tREF\n").encode()) # variant is only present in imputed files, non-zero genotype [(1, 0), (1, 1), (0, 1)].
                                else:
                                    fw.write((f"{chrom}\t{imp_pos}\t{imp_ref}\t{imp_alt}\t{imp_gt}\t{None}\tREF_0ALT\n").encode()) # variant is only present in imputed files, zero genotype [(0,0)].
                        else:
                            break    
                    if (imputed_truth == False):
                        if(truth_gt != 0):
                            fw.write((f"{chrom}\t{truth_pos}\t{truth_ref}\t{truth_alt}\t{None}\t{truth_gt}\tWGS\n").encode()) # variant is only present in truth files, non-zero genotype [(1, 0), (1, 1), (0, 1)].
                        else:
                            fw.write((f"{chrom}\t{truth_pos}\t{truth_ref}\t{truth_alt}\t{None}\t{truth_gt}\tWGS_0ALT\n").encode()) # variant is only present in truth files, zero genotype [(0, 0)].

                for imp_pos, imp_ref, imp_alt, imp_gt in imp_variants:
                        imp_variants_buffer.append((imp_pos, imp_ref, imp_alt, imp_gt))      
                for imp_pos, imp_ref, imp_alt, imp_gt in imp_variants_buffer:
                        if(imp_gt != 0):
                            fw.write((f"{chrom}\t{imp_pos}\t{imp_ref}\t{imp_alt}\t{imp_gt}\t{None}\tREF\n").encode()) # variant is only present in imputed files, non-zero genotype [(1, 0), (1, 1), (0, 1)].
                        else:
                            fw.write((f"{chrom}\t{imp_pos}\t{imp_ref}\t{imp_alt}\t{imp_gt}\t{None}\tREF_0ALT\n").encode()) # variant is only present in imputed files, zero genotype [(0, 0)].

if __name__ == "__main__":
    args = argparser.parse_args()
    sample_name = args.in_sample_name   
    path_imputed = args.in_imp_vcf
    path_truth = args.in_truth_vcf
    path_out = args.out_file_name

    compare(path_imputed, path_truth, sample_name, path_out)
 
