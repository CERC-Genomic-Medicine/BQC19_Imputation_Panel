#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/

process get_ref_chr_names{
    cache "lenient"

    input:
    tuple val(chrX), val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple stdout, val(chrX), val(sex_id), path(vcf), path(vcf_index)

    """
    tabix -l ${vcf} 
    """
}

process rm_chr_name_ref{
    cache "lenient"
    
    input:
    tuple val(chr_name), val(chrX), val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple val(chr_name), val(chrX), val(sex_id), path("*.vcf.gz"), path("*.vcf.gz.tbi")
    
    if (chrX == False) {
    """
    echo "${chr_name}" > chroms1.txt
    chr=${chr_name}
    echo "\${chr:3}" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.vcf.gz
    """
    } else {
    """
    echo "${chr_name}" > chroms1.txt
    echo "\$X" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.vcf.gz
    """
    }
}

process convert_ref_vcf{
    cache "lenient"

    input:
    tuple val(chr_name), val(chrX), val(sex_id), path(ref_vcf), path(ref_vcf_index)
    
    output:
    tuple val(chr_name), val(chrX), val(sex_id), path( "*.m3vcf.gz")
    
    publishDir "minimac_m3vcfs/", pattern: "*.m3vcf.gz", mode: "copy"

    """
    ${params.minimac3} --refHaps ${ref_vcf.getBaseName()}.ref.vcf.gz  --processReference --prefix ${ref_vcf.getBaseName()}
    """
}

process get_study_chr_names{
    cache "lenient"

    input:
    tuple val(chrX), val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple stdout, val(chrX), val(sex_id), path(vcf), path(vcf_index)

    """
    tabix -l ${vcf}
    """
}

process rm_chr_name_study{
    cache "lenient"
    
    input:
    tuple val(chr_name), val(chrX), val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple val(chr_name), val(chrX), val(sex_id), path("*.vcf.gz"), path("*.vcf.gz.tbi")
    
    if (chrX == False) {
    """
    echo "${chr_name}" > chroms1.txt
    chr=${chr_name}
    echo "\${chr:3}" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.vcf.gz
    """
    } else {
    """
    echo "${chr_name}" > chroms1.txt
    echo "\$X" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.vcf.gz
    """
    }
}

process minimac_imputation{
    cache "lenient"
    
    input:
    tuple val(chr_name), val(chrX), val(sex_id), file(ref_vcf), file(study_vcf), file(study_vcf_index)
     
    output:
    tuple path("*.vcf.gz*"), path("*.info")

    publishDir "imputed_vcfs/", pattern: "*.vcf.gz*", mode: "copy"
    publishDir "imputed_info/", pattern: "*.info", mode: "copy"
    
    if (chrX == False) {
    """
    ${params.minimac4} --refHaps $ref_vcf --haps $study_vcf --prefix ${study_vcf.getBaseName()} --meta --ignoreDuplicates
    echo "${chr_name}" > chroms1.txt
    chr=${chr_name}
    echo "\${chr:3}" > chroms2.txt

    paste chroms2.txt chroms1.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${study_vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${study_vcf.getBaseName()}.vcf.gz
    """
    } else {
    """
    ${params.minimac4} --refHaps $ref_vcf --haps $study_vcf --prefix ${study_vcf.getBaseName()} --meta --ignoreDuplicates

    echo "X" > chroms1.txt
    echo "chrX" > chroms2.txt

    paste chroms2.txt chroms1.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${study_vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${study_vcf.getBaseName()}.vcf.gz
    """
    }
}
workflow {

        ref_ch = Channel.fromPath(params.ref_vcf_path).map{ vcf -> [ vcf.name.toString().tokenize('.').contains('chrX'), vcf.name.toString().tokenize('.').contains('female'), vcf, vcf + ".tbi" ] }
        study_ch = Channel.fromPath(params.study_vcf_path).map{ vcf -> [ vcf.name.toString().tokenize('.').contains('chrX'), vcf.name.toString().tokenize('.').contains('female'), vcf, vcf + ".tbi" ] }

        ref_vcfs = get_ref_chr_names(ref_ch)
        study_vcfs = get_study_chr_names(study_ch)
        
        ref_rm_chr_vcfs = rm_chr_name_ref(ref_vcfs)
        study_rm_chr_vcfs = rm_chr_name_study(study_vcfs)

        ref_cnv_vcfs = convert_ref_vcf(ref_rm_chr_vcfs)

        imputation_ch = ref_cnv_vcfs.join(study_rm_chr_vcfs, by:[0, 1, 2])

        minimac_imputation(imputation_ch)
}