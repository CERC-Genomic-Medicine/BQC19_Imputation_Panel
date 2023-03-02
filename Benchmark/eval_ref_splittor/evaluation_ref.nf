#!/usr/bin/env nextflow


input_ch = Channel.fromPath( params.vcf_path )



process extract_ref_chr_num{
    cache "lenient"

    input:
    set file(vcf), file(vcf_index) from input_ch.map{ vcf -> [ vcf, vcf + ".tbi" ] }

    output:    
    tuple file("*.ref.vcf.gz"), file("*.ref.vcf.gz.tbi") into ref_results mode flatten
    tuple file("*.evaluation.vcf.gz"), file("*.evaluation.vcf.gz.tbi") into evaluation_results mode flatten

    publishDir "ref_vcfs/", pattern: "*.ref.vcf.gz", mode: "copy"
    publishDir "ref_vcfs/", pattern: "*.ref.vcf.gz.tbi", mode: "copy"

    publishDir "evaluation_vcfs/", pattern: "*.evaluation.vcf.gz", mode: "copy"
    publishDir "evaluation_vcfs/", pattern: "*.evaluation.vcf.gz.tbi", mode: "copy"


    """
    bcftools view -S  ${params.unrelated} --force-samples $vcf -Oz -o ${vcf.getBaseName()}.evaluation.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.evaluation.vcf.gz
    bcftools view -S  ^${params.unrelated} --force-samples $vcf -Oz -o ${vcf.getBaseName()}.ref.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.ref.vcf.gz
    """
}

