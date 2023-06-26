/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>
* VERSION: 1.0
* YEAR: 2022
*/
process imputed_vs_imputed {
   //errorStrategy "retry"
   //maxRetries 3
   cache "lenient"
   cpus 1
   memory "128GB"
   time "3h"
   //scratch true

   input:
   tuple val(sample_name), val(first_reference_name), path(first_post_imputation_file), val(second_reference_name), path(second_post_imputation_file)

   output:
   path("*_analysis.txt"), emit: analysis_files 
   path("*_pm.txt"), emit: present_one_missing_other_files
   path("*_bw.txt"), emit: well_one_bad_other_files
   path("*_bm.txt"), emit: bad_one_missing_other_files
   path("*_wm.txt"), emit: well_one_missing_other_files
   path("*_ww.txt"), emit: well_one_well_other_files
   path("*_mm.txt"), emit: missing_one_missing_other_files
   path("*_bb.txt"), emit: bad_one_bad_other_files



   publishDir "new_result/post_imputation_analysis/", pattern: "*_analysis.txt", mode: "copy"
   publishDir "new_result/bad_one_well_other/", pattern: "*_bw.txt", mode: "copy"
   publishDir "new_result/present_one_missing_other/", pattern: "*_pm.txt", mode: "copy"
   publishDir "new_result/bad_one_missing_other/", pattern: "*_bm.txt", mode: "copy"   
   publishDir "new_result/well_one_missing_other/", pattern: "*_wm.txt", mode: "copy"
   publishDir "new_result/well_one_well_other/", pattern: "*_ww.txt", mode: "copy"
   publishDir "new_result/missing_one_missing_other/", pattern: "*_mm.txt", mode: "copy"
   publishDir "new_result/bad_one_bad_other/", pattern: "*_bb.txt", mode: "copy"

    """
    imputed_vs_imputed.py -fq ${first_post_imputation_file} -sq ${second_post_imputation_file}  -s ${sample_name}  -fr ${first_reference_name} -sr ${second_reference_name} 
    """
}


workflow {
    first_reference_file = Channel.fromPath(params.post_imputation_files_first_reference).map{ file -> [file.name.toString().tokenize('_').get(1), file.name.toString().tokenize('_').get(0), file] }

    second_reference_file = Channel.fromPath(params.post_imputation_files_second_reference).map{ file -> [file.name.toString().tokenize('_').get(1), file.name.toString().tokenize('_').get(0), file] }

    combine_channel = first_reference_file.join(second_reference_file)
    files = imputed_vs_imputed(combine_channel)    
}