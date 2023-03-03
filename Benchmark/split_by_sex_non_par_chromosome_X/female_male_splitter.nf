#!/usr/bin/env nextflow

chr_X = Channel.fromPath( params.chromosome_X_path )

process split_by_sex{

    cache "lenient"

    input:
    set file(vcf), file(vcf_index) from chr_X.map{ vcf -> [ vcf, vcf + ".tbi" ] }

    output:
    tuple file("*.vcf.gz"), file("*.vcf.gz.tbi") into out_ch mode flatten
    
    publishDir "splitted_X_female_male/", pattern: "*.vcf.gz", mode: "copy"
    publishDir "splitted_X_female_male/", pattern: "*.vcf.gz.tbi", mode: "copy"

    """
    bcftools view  $vcf -S ${params.female} --force-samples -Oz -o ${vcf.getBaseName()}.female.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.female.vcf.gz
    
    bcftools view  $vcf -S ${params.male} --force-samples -Oz -o ${vcf.getBaseName()}.male.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.male.vcf.gz
    """
}