/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>
* VERSION: 1.0
* YEAR: 2023
*/

process metaimputation{
    cache "lenient"
    
    input:
    set val(ref_chr_num), val(first_ref_name), val(second_ref_name)
     
    output:
    file("*.vcf.gz") 

    publishDir "meta_imputed_vcfs/", pattern: "*.vcf.gz", mode: "copy"
    """
    first_ref_file_name=${params.first_ref_path}${ref_chr_num}${first_ref_name}
    second_ref_file_name=${params.second_ref_path}${ref_chr_num}${second_ref_name}
    ${params.metaminimac2} -i \${first_ref_file_name}:\${second_ref_file_name} --skipPhasingCheck "[ON]"
    """
}

workflow {
   first_imputed_files = ref_ch = Channel.fromPath(params.first_ref_files).map { file -> [ file.name.toString().tokenize('.').get(1), file.name.toString().tokenize('.').get(0)] }

   second_imputed_files = Channel.fromPath(params.second_ref_files).map { file -> [ file.name.toString().tokenize('.').get(1), file.name.toString().tokenize('.').get(0)] }

   join_by_chr = first_imputed_files.join(second_imputed_files)

   metaimputation(join_by_chr)
   
 }