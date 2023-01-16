/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>
* VERSION: 1.0
* YEAR: 2022
*/

process get_imputed_chr_names {
   cache "lenient"
   cpus 1
   memory "4GB"
   time "1h"
   scratch true

   input:
   tuple path(vcf), path(vcf_index)

   output:
   tuple stdout, path(vcf), path(vcf_index)

   """
   tabix -l ${vcf}
   """
} 

process get_truth_chr_names {
   cache "lenient"
   cpus 1
   memory "4GB"
   time "1h"
   scratch true

   input:
   tuple path(vcf), path(vcf_index)

   output:
   tuple stdout, path(vcf), path(vcf_index)

   """
   tabix -l ${vcf}
   """
}

process get_imputed_sample_names {
   cache "lenient"
   cpus 1
   memory "4GB"
   time "1h"
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
   cache "lenient"
   cpus 1
   memory "4GB"
   time "1h"
   scratch true

   input:
   tuple val(chromosome), path(imputed_vcf), path(imputed_vcf_index), path(truth_vcf), path(truth_vcf_index)
   each individual

   output:
   path("*.txt")

   publishDir "result/${params.ref_name}/", pattern: "*_post_imputation_analysis.gz", mode: "copy"

    """
    imputed_vs_truth2.py -iv ${imputed_vcf} -tv ${truth_vcf} -s ${individual} -o "${params.ref_name}_${individual}_${chromosome}_post_imputation_analysis.gz"
    """
}


workflow {
    imputed_files = Channel.fromPath(params.imputed_files).map{ vcf -> [ vcf, vcf + ".tbi" ] }
    truth_files = Channel.fromPath(params.truth_files).map{ vcf -> [ vcf, vcf + ".tbi" ] }
    
    imputed_by_chr = get_imputed_chr_names(imputed_files).map{ it -> [it[0].startsWith("chr") ? it[0].substring(3).trim() : it[0].trim(), it[1], it[2]] }
    truth_by_chr = get_truth_chr_names(truth_files).map{ it -> [it[0].startsWith("chr") ? it[0].substring(3).trim() : it[0].trim(), it[1], it[2]] }

    all_by_chr = imputed_by_chr.join(truth_by_chr)
    imputed_sample_names = get_imputed_sample_names(Channel.fromPath(params.imputed_files).first().map{ vcf -> [vcf, vcf + ".tbi"] }).flatMap{ it -> it.split("\t") }

    //all_by_chr.view()
    
    imputed_vs_truth(all_by_chr, imputed_sample_names)
}
