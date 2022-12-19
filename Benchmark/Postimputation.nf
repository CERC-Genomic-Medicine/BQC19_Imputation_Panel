/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>
* VERSION: 1.0
* YEAR: 2022
*/
process extract_chr_num {
   errorStrategy "retry"
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16GB"
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

process extract_chr_num2 {
   errorStrategy "retry"
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16GB"
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

process imputed_vs_truth {
   errorStrategy "retry"
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16GB"
   time "1h"
   scratch true

   input:
   tuple val(chromosome), path(imputed_vcf), path(imputed_vcf_index), path(truth_vcf), path(truth_vcf_index)
   each individual

   output:

    publishDir "result/${individual}", pattern: "*.txt", mode: "copy"

    """
    ImputedVsTruth.py -iv ${imputed_vcf} -tv ${truth_vcf} -s ${individual} -r ${ref_name} -c ${chromosome} -o ${individual}_${ref_name}_${chromosome}
    """
}


workflow {
    imputed_files = Channel.fromPath(params.imputed_files_first_reference).map{ vcf -> [ vcf, vcf + ".tbi" ] }
    truth_files = Channel.fromPath(params.truth_files).map{ vcf -> [ vcf, vcf + ".tbi" ] }
    samples = Channel.from(file(params.test_files).readLines())
    truth_chromosome_files = extract_chr_num(truth_files)
    imputed_chromosome_files = extract_chr_num2(imputed_files)
    imputed_truth_chromosome_files = imputed_chromosome_files.join(truth_chromosome_files)
    imputed_vs_truth(imputed_truth_chromosome_files, samples)
}