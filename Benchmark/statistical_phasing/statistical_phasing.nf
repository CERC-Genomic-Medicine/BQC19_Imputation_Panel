#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/


process beagle_statistical_phasing {
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    cpus 8
    memory "32GB"
    time "5h"
    input:
    tuple path(chr), path(index)
    
    output:
    tuple path("*.ref.vcf.gz"), path("*.ref.log")
    
    publishDir "phased/ref_vcfs/", pattern: "*.vcf.gz", mode: "copy"
    publishDir "phased/ref_logs/", pattern: "*.ref.log", mode: "copy"
    
    """
    java -jar -Djava.io.tmpdir=./temp/ -Xmx32g ${params.beagle} window=25.0 overlap=2.5 nthreads=8 gt=$chr out=${chr.getBaseName()}.ref 
    """
}

process remove_singletons {
    cache "lenient"
    cpus 1
    memory "16GB"
    time "5h"
    input:
    tuple path(chr), path(log)
    
    output:
    tuple path("*.vcf.gz"), path("*.vcf.gz.tbi")
    
    publishDir "phased/ref_withoutsingletons_vcfs/", pattern: "*.vcf.gz", mode: "copy"
    publishDir "phased/ref_withoutsingletons_vcfs/", pattern: "*.vcf.gz.tbi", mode: "copy"

    """
    bcftools view $chr -c 2 -Oz -o  ${chr.getBaseName()}.with_out_singletons.vcf.gz
    bcftools index --tbi ${chr.getBaseName()}.with_out_singletons.vcf.gz
    """
}

workflow {

        stat_phasing_ch = Channel.fromPath(params.ref_vcf_path).map{ vcf -> [vcf, vcf + ".tbi" ] }

        phased_vcfs = beagle_statistical_phasing(stat_phasing_ch)

        remove_singletons(phased_vcfs)


}
