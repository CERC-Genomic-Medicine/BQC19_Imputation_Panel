/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>
* VERSION: 1.0
* YEAR: 2023
*/

process metaimputation{
    cache "lenient"
    
    input:
    tuple val(ref_chr_num), val(first_ref_name), path(first_empirical_dose), path(first_dose), val(second_ref_name), path(second_empirical_dose), path(second_dose)
     
    output:
    file("*.gz") 

    publishDir "meta_imputed_vcfs/", pattern: "*.gz", mode: "copy"
    """
    first_ref_file_name=./${first_ref_name}.${ref_chr_num}
    second_ref_file_name=./${second_ref_name}.${ref_chr_num}
    ${params.metaminimac2} -i \${first_ref_file_name}:\${second_ref_file_name} --skipPhasingCheck "[ON]" --weight "[ON]" -o Meta.${first_ref_name}.${second_ref_name}.${ref_chr_num}
    """
}

workflow {
   first_imputed_files = Channel.fromPath(params.first_ref_files).map { file -> [ file.name.toString().tokenize('.').get(1), file.name.toString().tokenize('.').get(0), params.first_ref_path + file.name.toString().tokenize('.').get(0)+'.'+ file.name.toString().tokenize('.').get(1)+'.empiricalDose.vcf.gz', params.first_ref_path + file.name.toString().tokenize('.').get(0)+'.'+file.name.toString().tokenize('.').get(1)+'.dose.vcf.gz'] }

   second_imputed_files = Channel.fromPath(params.second_ref_files).map { file -> [ file.name.toString().tokenize('.').get(1), file.name.toString().tokenize('.').get(0), params.second_ref_path + file.name.toString().tokenize('.').get(0)+'.'+ file.name.toString().tokenize('.').get(1)+'.empiricalDose.vcf.gz', params.second_ref_path + file.name.toString().tokenize('.').get(0)+'.'+file.name.toString().tokenize('.').get(1)+'.dose.vcf.gz'] }

   join_by_chr = first_imputed_files.join(second_imputed_files)

   metaimputation(join_by_chr)
   
}