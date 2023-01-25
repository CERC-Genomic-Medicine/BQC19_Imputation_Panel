#!/usr/bin/env nextflow


chr_ch = Channel.fromPath(params.vcfs)


process beagle_statistical_phasing {

    input:
    file chr from chr_ch
    
    output:
    tuple file ("*.ref.vcf.gz"), file("*.ref.log") into ref_ch mode flatten
    
    publishDir "phased/ref_vcfs/", pattern: "*.vcf.gz", mode: "copy"
    publishDir "phased/ref_logs/", pattern: "*.ref.log", mode: "copy"
    
    """
    java -jar -Djava.io.tmpdir=./temp/ -Xmx32g ${params.beagle} window=25.0 overlap=2.5 nthreads=8 gt=$chr out=${chr.getBaseName()}.ref 
    """
}

process remove_singletons {
    input:
    set file(chr), file(log) from ref_ch
    
    output:
    tuple file("*.vcf.gz"), file("*.vcf.gz.tbi") into rm_singletons_ch
    
    publishDir "phased/ref_withoutsingletons_vcfs/", pattern: "*.vcf.gz", mode: "copy"
    publishDir "phased/ref_withoutsingletons_vcfs/", pattern: "*.vcf.gz.tbi", mode: "copy"

    """
    bcftools view $chr -c 2 -Oz -o  ${chr.getBaseName()}.with_out_singletons.vcf.gz
    bcftools index --tbi ${chr.getBaseName()}.with_out_singletons.vcf.gz
    """
}
