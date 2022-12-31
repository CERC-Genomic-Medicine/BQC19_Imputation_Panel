/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>
* VERSION: 1.0
* YEAR: 2022
*/
process imputed_vs_imputed {
   errorStrategy "retry"
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16GB"
   time "1h"
   scratch true

   input:
   tuple val(sample_name), val(chromosome), val(first_reference_name), path(first_post_imputation_file), val(second_reference_name), path(second_post_impputation_file)

   output:
   path("*.txt")

   publishDir "result/${params.mode}/post_imputation_analysis/", pattern: "*_post_imputation_analysis.txt", mode: "copy"
   publishDir "result/${params.mode}/bad_one_well_other/", pattern: "*_bw.txt", mode: "copy"
   publishDir "result/${params.mode}/well_one_missing_other/", pattern: "*_wm.txt", mode: "copy"
   publishDir "result/${params.mode}/missing_both/", pattern: "*_mm.txt", mode: "copy"

    """
    imputed_vs_imputed.py -fq ${first_post_imputation_file} -sq ${second_post_imputation_file}  -s ${sample_name}  -fr ${first_reference_name} -sr ${second_reference_name} -c ${chromosome} -m ${params.mode}
    """
}


workflow {
    first_reference_file = Channel.fromPath(params.post_imputation_files_first_reference).\
    map{ file -> file.name.toString().tokenize('_')[0:3] + [file ] }

    second_reference_file = Channel.fromPath(params.post_imputation_files_second_reference).\
    map{ file -> file.name.toString().tokenize('_')[0:3] + [file ] }

    combine_channel = first_reference_file.join(second_reference_file, by = [0, 1])

    
}