#!/usr/bin/env nextflow


raw_ch = Channel.fromPath( params.vcf_path )
panel_ch = Channel.fromPath( params.panel_path )



process extract_ref_chr_num{
    cache "lenient"

    input:
    set file(vcf), file(vcf_index) from raw_ch.map{ vcf -> [ vcf, vcf + ".tbi" ] }
    set file(vcf2), file(vcf2_index) from panel_ch.map{ vcf -> [ vcf, vcf + ".tbi" ] }

    output:    
    tuple file("*.ref.vcf.gz"), file("*.ref.vcf.gz.tbi") into ref_results mode flatten
    tuple file("*.evaluation.vcf.gz"), file("*.evaluation.vcf.gz.tbi") into evaluation_results mode flatten

    publishDir "ref_vcfs/ex1/", pattern: "*.ex1.ref.vcf.gz", mode: "copy"
    publishDir "ref_vcfs/ex1/", pattern: "*.ex1.ref.vcf.gz.tbi", mode: "copy"

    publishDir "evaluation_vcfs/ex1/", pattern: "*.ex1.evaluation.vcf.gz", mode: "copy"
    publishDir "evaluation_vcfs/ex1/", pattern: "*.ex1.evaluation.vcf.gz.tbi", mode: "copy"

    publishDir "ref_vcfs/ex2/", pattern: "*.ex2.ref.vcf.gz", mode: "copy"
    publishDir "ref_vcfs/ex2/", pattern: "*.ex2.ref.vcf.gz.tbi", mode: "copy"

    publishDir "evaluation_vcfs/ex2/", pattern: "*.ex2.evaluation.vcf.gz", mode: "copy"
    publishDir "evaluation_vcfs/ex2/", pattern: "*.ex2.evaluation.vcf.gz.tbi", mode: "copy"

    publishDir "ref_vcfs/ex3/", pattern: "*.ex3.ref.vcf.gz", mode: "copy"
    publishDir "ref_vcfs/ex3/", pattern: "*.ex3.ref.vcf.gz.tbi", mode: "copy"

    publishDir "evaluation_vcfs/ex3/", pattern: "*.ex3.evaluation.vcf.gz", mode: "copy"
    publishDir "evaluation_vcfs/ex3/", pattern: "*.ex3.evaluation.vcf.gz.tbi", mode: "copy"


    publishDir "ref_vcfs/ex4/", pattern: "*.ex4.ref.vcf.gz", mode: "copy"
    publishDir "ref_vcfs/ex4/", pattern: "*.ex4.ref.vcf.gz.tbi", mode: "copy"

    publishDir "evaluation_vcfs/ex4/", pattern: "*.ex4.evaluation.vcf.gz", mode: "copy"
    publishDir "evaluation_vcfs/ex4/", pattern: "*.ex4.evaluation.vcf.gz.tbi", mode: "copy"


    publishDir "ref_vcfs/ex5/", pattern: "*.ex5.ref.vcf.gz", mode: "copy"
    publishDir "ref_vcfs/ex5/", pattern: "*.ex5.ref.vcf.gz.tbi", mode: "copy"

    publishDir "evaluation_vcfs/ex5/", pattern: "*.ex5.evaluation.vcf.gz", mode: "copy"
    publishDir "evaluation_vcfs/ex5/", pattern: "*.ex5.evaluation.vcf.gz.tbi", mode: "copy"

    """
    bcftools view -S  ${params.ex1} --force-samples $vcf -Oz -o ${vcf.getBaseName()}.ex1.evaluation.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.ex1.evaluation.vcf.gz
    bcftools view -S  ^${params.ex1} --force-samples $vcf2 -Oz -o ${vcf2.getBaseName()}.ex1.ref.vcf.gz
    bcftools index --tbi ${vcf2.getBaseName()}.ex1.ref.vcf.gz

    bcftools view -S  ${params.ex2} --force-samples $vcf -Oz -o ${vcf.getBaseName()}.ex2.evaluation.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.ex2.evaluation.vcf.gz
    bcftools view -S  ^${params.ex2} --force-samples $vcf2 -Oz -o ${vcf2.getBaseName()}.ex2.ref.vcf.gz
    bcftools index --tbi ${vcf2.getBaseName()}.ex2.ref.vcf.gz

    bcftools view -S  ${params.ex3} --force-samples $vcf -Oz -o ${vcf.getBaseName()}.ex3.evaluation.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.ex3.evaluation.vcf.gz
    bcftools view -S  ^${params.ex3} --force-samples $vcf2 -Oz -o ${vcf2.getBaseName()}.ex3.ref.vcf.gz
    bcftools index --tbi ${vcf2.getBaseName()}.ex3.ref.vcf.gz

    bcftools view -S  ${params.ex4} --force-samples $vcf -Oz -o ${vcf.getBaseName()}.ex4.evaluation.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.ex4.evaluation.vcf.gz
    bcftools view -S  ^${params.ex4} --force-samples $vcf2 -Oz -o ${vcf2.getBaseName()}.ex4.ref.vcf.gz
    bcftools index --tbi ${vcf2.getBaseName()}.ex4.ref.vcf.gz

    bcftools view -S  ${params.ex5} --force-samples $vcf -Oz -o ${vcf.getBaseName()}.ex5.evaluation.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.ex5.evaluation.vcf.gz
    bcftools view -S  ^${params.ex5} --force-samples $vcf2 -Oz -o ${vcf2.getBaseName()}.ex5.ref.vcf.gz
    bcftools index --tbi ${vcf2.getBaseName()}.ex5.ref.vcf.gz

    """
}

