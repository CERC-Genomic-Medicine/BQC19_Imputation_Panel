#!/usr/bin/env nextflow


chr_ch = Channel.fromPath(params.vcfs)


process preprocessing {

    input:
    file chr from chr_ch
    
    output:
    tuple file ("*.vcf.gz"), file("*.vcf.gz.tbi") into ref_ch
    
    publishDir "preprocessed_vcfs/", pattern: "*.vcf.gz", mode: "copy"   
    publishDir "preprocessed_vcfs/", pattern: "*.vcf.gz.tbi", mode: "copy"    

    """
    bcftools view -f PASS -m2 -M2 -v snps,indels -i'INFO/AN>4183' $chr | bcftools annotate -x ^INFO/AF,^INFO/AN,^INFO/AC,^FORMAT/GT | bcftools norm -d all -Oz -o ${chr.getBaseName()}.filtered.vcf.gz
    bcftools tabix --tbi ${chr.getBaseName()}.filtered.vcf.gz
    """
}
