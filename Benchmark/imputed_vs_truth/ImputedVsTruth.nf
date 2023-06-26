/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>
* VERSION: 1.0
* YEAR: 2023
*/

process concat_imputed {
   cache "lenient"
   cpus 1
   memory "64GB"
   time "04:00:00"
   scratch true
   input:
   path(vcfs)
   output:
   tuple path("*.vcf.gz"), path("*.vcf.gz.tbi")
   publishDir "result/${params.ref_name}/concat/", pattern: "*.gz", mode: "copy"
   publishDir "result/${params.ref_name}/concat/", pattern: "*.gz.tbi", mode: "copy"

   """
   bcftools concat ${vcfs} -Oz -o imputed.vcf.gz
   bcftools index --tbi imputed.vcf.gz
   """
} 
process concat_truth {
   cache "lenient"
   cpus 1
   memory "64GB"
   time "04:00:00"
   scratch true
   input:
   path(vcfs)
   output:
   tuple path("*.vcf.gz"), path("*.vcf.gz.tbi")
   publishDir "result/${params.ref_name}/concat/", pattern: "*.gz", mode: "copy"
   publishDir "result/${params.ref_name}/concat/", pattern: "*.gz.tbi", mode: "copy"

   """
   bcftools concat ${vcfs} -Oz -o truth.vcf.gz
   bcftools index --tbi truth.vcf.gz
   """
}
process get_imputed_sample_names {
   cache "lenient"
   cpus 1
   memory "4GB"
   time "00:30:00"
   scratch true
   input:
   tuple path(vcf), path(vcf_index)
   output:
   stdout
   """
   tabix -H ${vcf} | tail -n1 | cut -f10-
   """
}
process imputed_vs_truth {
   //errorStrategy "retry"
   //maxRetries 3
   //debug true
   cache "lenient"
   cpus 1
   memory "16GB"
   time "07:00:00"
   scratch true
   input:
   tuple path(imputed_vcf), path(imputed_vcf_index)
   tuple path(truth_vcf), path(truth_vcf_index)
    each individual

    output:
    tuple val(individual), path("*${individual}*.txt.gz")

    publishDir "result/${params.ref_name}/all/", pattern: "*.txt.gz", mode: "copy"

     """
     imputed_vs_truth.py -iv ${imputed_vcf} -tv ${truth_vcf} -s ${individual} -o "${params.ref_name}_${individual}_post_imputation_analysis.txt.gz"
     """
 }

 
process generate_summary {
   cache "lenient"
   cpus 1
   memory "16GB"
   time "00:30:00"
   scratch true

   input:
   tuple val(individual), path(quality_file_per_sample)

   output:
   path("*.txt")

   publishDir "result/${params.ref_name}/summary/", pattern: "*.txt", mode: "copy"
   """
   generate_summary.py -i ${quality_file_per_sample} -s ${individual} -o ${individual}.summary.txt
   """
}

process concat_all_samples_summary {
   cache "lenient"
   cpus 1
   memory "4GB"
   time "00:15:00"
   scratch true

   input:
   path(summary_per_sample)

   output:
   path("*.txt")
   publishDir "result/${params.ref_name}/summary_concat/", pattern: "*.txt", mode: "copy"

   """
   awk 'FNR>1' ${summary_per_sample} > ${params.ref_name}_concat_all_summary.txt
   sed -i -e '1iSample ID\tWGS\tREF\tWGS_AND_REF\tWGS_AND_REF_EQ\tWGS_AND_REF_LT\tWGS_AND_REF_GT\tREF_0ALT\tWGS_0ALT\t'AA_Concordance_FRAC'\t'AA_Concordance'\t'Coverage'\t'RA_Discordance_FRAC'\t'RA_Discordance'\n' ${params.ref_name}_concat_all_summary.txt
   """
}


workflow {
   imputed_files = Channel.fromPath(params.imputed_files)
   truth_files = Channel.fromPath(params.truth_files)

   imputed_sample_names = get_imputed_sample_names(Channel.fromPath(params.truth_files).first().map{ vcf -> [vcf, vcf + ".tbi"] }).flatMap{ it -> it.split("\n")}.flatMap{ it -> it.split("\t")}
   imputed_vcf = concat_imputed(imputed_files.collect())
   truth_vcf = concat_truth(truth_files.collect())
   quality_files = imputed_vs_truth(imputed_vcf, truth_vcf, imputed_sample_names)

   summary_files = generate_summary(quality_files) 
   concat_all_samples_summary(summary_files.collect())

 }