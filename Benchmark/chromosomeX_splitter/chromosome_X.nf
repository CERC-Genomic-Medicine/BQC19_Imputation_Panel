#!/usr/bin/env nextflow

chr_X = Channel.fromPath( params.chromosome_X_path )

process extract_PAR_regions_chr_X_annotate{

    cache "lenient"

    input:
    set file(vcf), file(vcf_index) from chr_X.map{ vcf -> [ vcf, vcf + ".tbi" ] }

    output:
    tuple file("*.PAR1_23.vcf.gz"), file("*.PAR1_23.vcf.gz.tbi") into PAR1 mode flatten
    tuple file("*.Non_PAR_24.vcf.gz"), file("*.Non_PAR_24.vcf.gz.tbi") into non_PAR mode flatten
    tuple file("*.PAR2_25.vcf.gz"), file("*.PAR2_25.vcf.gz.tbi") into PAR2 mode flatten
    publishDir "splitted_X/", pattern: "*.vcf.gz", mode: "copy"
    publishDir "splitted_X/", pattern: "*.vcf.gz.tbi", mode: "copy"

    """
    bcftools view  $vcf -r chrX:1-2781479 -Oz -o PAR1.vcf.gz
    bcftools annotate --rename-chrs $params.PAR1_rename PAR1.vcf.gz -Oz -o ${vcf.getBaseName()}.PAR1_23.vcf.gz 
    bcftools index --tbi ${vcf.getBaseName()}.PAR1_23.vcf.gz
    bcftools view  $vcf -r chrX:2781479-155701383 -Oz -o Non_PAR.vcf.gz
    bcftools annotate --rename-chrs $params.Non_PAR_rename Non_PAR.vcf.gz -Oz -o ${vcf.getBaseName()}.Non_PAR_24.vcf.gz 
    bcftools index --tbi ${vcf.getBaseName()}.Non_PAR_24.vcf.gz
    bcftools view  $vcf -r chrX:155701383-156040896 -Oz -o PAR2.vcf.gz
    bcftools annotate --rename-chrs $params.PAR2_rename PAR2.vcf.gz -Oz -o ${vcf.getBaseName()}.PAR2_25.vcf.gz 
    bcftools index --tbi ${vcf.getBaseName()}.PAR2_25.vcf.gz
    """
}