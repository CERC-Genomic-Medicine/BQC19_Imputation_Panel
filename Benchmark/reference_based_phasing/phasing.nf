#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/

process get_ref_chr_names{
    cache "lenient"
    cpus 1
    memory "4GB"
    time "00:30:00"
    scratch true
    input:
    tuple val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple stdout, val(sex_id), path(vcf), path(vcf_index)

    """
    tabix -l ${vcf} 
    """
}


process get_study_chr_names{
    cache "lenient"
    cpus 1
    memory "4GB"
    time "00:30:00"
    scratch true
    input:
    tuple val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple stdout, val(sex_id), path(vcf), path(vcf_index)

    """
    tabix -l ${vcf}
    """
}



process eagle_phasing{
    cache "lenient"
    cpus 1
    memory "16GB"
    time "5h"
    scratch true
    input:
    tuple val(chromosome), val(sex_id), file(ref_vcf), file(ref_vcf_index), file(study_vcf), file(study_vcf_index) 
        
    output:
    file ("*.phased.annotated.vcf.gz")
    publishDir "phased_vcfs/", pattern: "*.vcf.gz", mode: "copy"

    script:
    if(params.chromosomeX == true){
    """
    echo "${chromosome}" > chroms1.txt
    chr=${chromosome}
    echo "chrX" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt   
    
    bcftools annotate --rename-chrs chr_name_conv.txt $study_vcf -Oz -o ${study_vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${study_vcf.getBaseName()}.vcf.gz

    ${params.eagle} --vcfRef $ref_vcf --vcfTarget ${study_vcf.getBaseName()}.vcf.gz --geneticMapFile ${params.genetic_map} --outPrefix ${study_vcf.getBaseName()}.phased --allowRefAltSwap --vcfOutFormat z

    paste chroms2.txt chroms1.txt > chr_name_conv.txt   

    bcftools annotate --rename-chrs chr_name_conv.txt ${study_vcf.getBaseName()}.phased.vcf.gz -Oz -o ${study_vcf.getBaseName()}.phased.annotated.vcf.gz
    bcftools index --tbi ${study_vcf.getBaseName()}.phased.annotated.vcf.gz

    """
    } else {
    """
    ${params.eagle} --vcfRef $ref_vcf --vcfTarget $study_vcf --geneticMapFile ${params.genetic_map} --outPrefix ${study_vcf.getBaseName()} --allowRefAltSwap --vcfOutFormat z
    """
    }

}

workflow {

        ref_ch = Channel.fromPath(params.ref_vcf_path).map{ vcf -> [ vcf.name.toString().tokenize('.').contains('female'), vcf, vcf + ".tbi" ] }
        study_ch = Channel.fromPath(params.study_vcf_path).map{ vcf -> [ vcf.name.toString().tokenize('.').contains('female'), vcf, vcf + ".tbi" ] }

        ref_vcfs = get_ref_chr_names(ref_ch)
        study_vcfs = get_study_chr_names(study_ch)

        phasing_ch = ref_vcfs.join(study_vcfs, by:[0, 1])

        eagle_phasing(phasing_ch)

}