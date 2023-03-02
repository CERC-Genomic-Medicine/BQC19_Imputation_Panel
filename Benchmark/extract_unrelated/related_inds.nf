#!/usr/bin/env nextflow


chr_ch = Channel.fromPath(params.vcfs)


process concat {

    input:
    file chr from chr_ch.collect()
    
    output:
    tuple file ("*.vcf.gz"), file("*.vcf.gz.tbi") into concat_ch
    
    publishDir "concat_vcf/", pattern: "*.vcf.gz", mode: "copy"   
    publishDir "concat_vcf/", pattern: "*.vcf.gz.tbi", mode: "copy"    

    """
    bcftools concat $chr -Oz -o all_chromosomes.vcf.gz
    bcftools tabix --tbi all_chromosomes.vcf.gz
    """
}

process remove_related{
    
     input:
     set file(all), file(all_index) from concat_ch
     
     output:
     tuple file("*.kin0"), file("*.log") into plink_ch
     
     publishDir "results/", pattern:"*.kin0", mode: "copy"     
     publishDir "results/", pattern:"*.log", mode: "copy"     
     """
     ${params.plink} --make-king-table --vcf  $all --king-table-filter 0.0442
     """
}

