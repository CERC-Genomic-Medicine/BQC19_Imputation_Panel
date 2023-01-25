#!/usr/bin/env nextflow


ref_ch = Channel.fromPath( params.ref_vcf_path )

study_ch = Channel.fromPath( params.study_vcf_path )

process extract_ref_chr_num{
    cache "lenient"

    input:
    set file(vcf), file(vcf_index) from ref_ch.map{ vcf -> [ vcf, vcf + ".tbi" ] }

    output:
    tuple stdout, file(vcf), file(vcf_index) into ref_vcfs mode flatten

    """
    tabix -l ${vcf} 
    """
}


process extract_study_chr_num{
    cache "lenient"

    input:
    set file(vcf), file(vcf_index) from study_ch.map{ vcf -> [ vcf, vcf + ".tbi" ] }

    output:
    set stdout, file(vcf), file(vcf_index) into study_vcfs mode flatten

    """
    tabix -l ${vcf}
    """
}


phasing_ch = ref_vcfs.join(study_vcfs)

process eagle_phasing{
    cache "lenient"
    
    input:
    set val(ref_chr_num), file(ref_vcf), file(ref_vcf_index), file(study_vcf), file(study_vcf_index) from phasing_ch
     
    output:
    file "*.vcf.gz" into phased_vcfs mode flatten
    publishDir "phased_vcfs/", pattern: "*.vcf.gz", mode: "copy"

    """
    ${params.eagle} --vcfRef $ref_vcf --vcfTarget $study_vcf --geneticMapFile ${params.genetic_map} --outPrefix ${study_vcf.getBaseName()} --allowRefAltSwap --vcfOutFormat z
    """
}