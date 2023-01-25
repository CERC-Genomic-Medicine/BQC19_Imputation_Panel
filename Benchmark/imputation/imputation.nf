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

process rm_chr_name_ref{
    cache "lenient"
    
    input:
    set val(chr_name), file(vcf), file(vcf_index) from ref_vcfs

    output:
    tuple val(chr_name), file("*.vcf.gz"), file("*.vcf.gz.tbi") into ref_comp_ch  mode flatten
    
    """
    echo "${chr_name}" > chroms1.txt
    echo "${chr_name:3}" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt   
    
    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.vcf.gz
    """
}



process convert_ref_vcf{
    cache "lenient"

    input:
    set val(ref_chr_num), file(ref_vcf), file(ref_vcf_index) from ref_comp_ch
    
    output:
    tuple val(ref_chr_num), file("*.m3vcf.gz") into converted_vcfs_ch mode flatten
    
    publishDir "minimac_m3vcfs/", pattern: "*.m3vcf.gz", mode: "copy"

    
    """
    ${params.minimac3} --refHaps $ref_vcf  --processReference --prefix ${ref_vcf.getBaseName()}
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

process rm_chr_name_study{
    cache "lenient"
    
    input:
    set val(chr_name), file(vcf), file(vcf_index) from study_vcfs

    output:
    tuple val(chr_name), file("*.vcf.gz"), file("*.vcf.gz.tbi") into study_comp_ch  mode flatten
    
    """
    echo "${chr_name}" > chroms1.txt
    echo "${chr_name:3}" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.vcf.gz
    """
}
imputation_ch = converted_vcfs_ch.join(study_comp_ch)

process imputation{
    cache "lenient"
    
    input:
    set val(ref_chr_num), file(ref_vcf), file(study_vcf), file(study_vcf_index) from imputation_ch
     
    output:
    tuple file("*.vcf.gz"), file("*.info") into imputed_results mode flatten

    publishDir "imputed_vcfs/", pattern: "*.vcf.gz", mode: "copy"
    publishDir "imputed_info/", pattern: "*.info", mode: "copy"
    """
    ${params.minimac4} --refHaps $ref_vcf --haps $study_vcf --prefix ${study_vcf.getBaseName()} --ignoreDuplicates
    """
}
