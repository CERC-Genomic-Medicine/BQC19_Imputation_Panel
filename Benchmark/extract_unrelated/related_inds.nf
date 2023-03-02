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
     tuple file("*.king.cutoff.in.id"), file("*.king.cutoff.out.id") into plink_ch
     
     publishDir "related_ind/", pattern:"*.king.cutoff.out.id", mode: "copy"
     publishDir "unrelated_ind/", pattern:"*.king.cutoff.in.id", mode: "copy"
     
     """
     ${params.plink} --make-king --vcf $all --king-cutoff 0.01105
     """


}