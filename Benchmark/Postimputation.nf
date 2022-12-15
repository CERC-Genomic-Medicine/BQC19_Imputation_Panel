/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>
* VERSION: 1.0
* YEAR: 2022
*/

process compute_depth {
   errorStrategy "retry"
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "4GB"
   time "1h"
   scratch true

   input:
   tuple path(bam), path(bam_index)
   each chromosome

   output:
   tuple val(chromosome),
      path("${chromosome}.${bam.getSimpleName()}.depth.gz"), 
      path("${chromosome}.${bam.getSimpleName()}.depth.gz.tbi")

    publishDir "result/depth/${chromosome}", pattern: "*.depth*", mode: "copy"

    """
    samtools depth -a -s -q20 -Q20 -r ${chromosome} ${bam} | bgzip > ${chromosome}.${bam.getSimpleName()}.depth.gz
    tabix -s1 -b2 -e2 ${chromosome}.${bam.getSimpleName()}.depth.gz
    """
}


process aggregate {
   errorStrategy "finish"
   cache "lenient"
   cpus 1
   memory "8GB"
   time "5d"
   scratch true

   input:
   tuple val(chromosome), path(depth_files), path(depth_indices)

   output:
   path("${chromosome}.aggregated.txt.gz")

   publishDir "result/aggregated/", pattern: "*.aggregated.txt.gz*", mode: "copy"

   """
   find . -name "${chromosome}.*.depth.gz" > files_list.txt
   aggregate.py -f t -i files_list.txt -o ${chromosome}.aggregated.txt.gz
   """
}


process summarize {
   errorStrategy "finish"
   cache "lenient"
   cpus 1
   memory "8GB"
   time "3h"
   scratch true

   input:
   path aggregate_file

   output:
   tuple path("${aggregate_file.getBaseName()}.PCT_INDV_OVER_${params.min_dp}X.summary.txt"), path(aggregate_file)
   
   publishDir "result/summary/", pattern: "*.txt", mode: "copy"

   """
   summarize.py -i ${aggregate_file} -c PCT_INDV_OVER_${params.min_dp}X -o ${aggregate_file.getBaseName()}.PCT_INDV_OVER_${params.min_dp}X.summary.txt
   """
}


process create_accessibility_mask {

   errorStrategy "finish"
   cache "lenient"
   cpus 1
   memory "4GB"
   time "1h"
   //scratch true

   input:
   tuple path(stats), path(aggregate_file)

   output:
   path("*.bed")

   publishDir "result/accessibility_mask/", pattern: "*.bed", mode: "copy"

   script:
   if ((params.min_pct_ind != null) && (params.max_mean_dp != null)) {
      """
      min_pct_ind=${params.min_pct_ind}   
      max_mean_dp=${params.max_mean_dp}
      accessibility_mask.py -i ${aggregate_file} -c PCT_INDV_OVER_${params.min_dp}X -m "\${min_pct_ind}" -M "\${max_mean_dp}" -o ${aggregate_file.getBaseName()}.PCT_INDV_OVER_${params.min_dp}X_\${min_pct_ind}_MEAN_\${max_mean_dp}.bed
      """
   } else if (params.min_pct_ind != null) {
      """
      min_pct_ind=${params.min_pct_ind}
      max_mean_dp=\$(grep "^99%" ${stats} | cut -f2)
      accessibility_mask.py -i ${aggregate_file} -c PCT_INDV_OVER_${params.min_dp}X -m "\${min_pct_ind}" -M "\${max_mean_dp}" -o ${aggregate_file.getBaseName()}.PCT_INDV_OVER_${params.min_dp}X_\${min_pct_ind}_MEAN_\${max_mean_dp}.bed
      """
   } else if (params.max_mean_dp != null) {
      """
      min_pct_ind=\$(grep "^5%" ${stats} | cut -f3)
      max_mean_dp=${params.max_mean_dp}
      accessibility_mask.py -i ${aggregate_file} -c PCT_INDV_OVER_${params.min_dp}X -m "\${min_pct_ind}" -M "\${max_mean_dp}" -o ${aggregate_file.getBaseName()}.PCT_INDV_OVER_${params.min_dp}X_\${min_pct_ind}_MEAN_\${max_mean_dp}.bed
      """
   } else {
      """
      min_pct_ind=\$(grep "^5%" ${stats} | cut -f3)
      max_mean_dp=\$(grep "^99%" ${stats} | cut -f2)
      accessibility_mask.py -i ${aggregate_file} -c PCT_INDV_OVER_${params.min_dp}X -m "\${min_pct_ind}" -M "\${max_mean_dp}" -o ${aggregate_file.getBaseName()}.PCT_INDV_OVER_${params.min_dp}X_\${min_pct_ind}_MEAN_\${max_mean_dp}.bed
      """
   }
}


workflow {
   if (params.compute_depth == true) {
       bam_files = Channel.fromPath(params.input_files).map{ file -> [file, file + (file.getExtension() == "bam" ? ".bai" : ".crai")] }
       chromosomes = Channel.from(params.chromosomes)
       depth_files = compute_depth(bam_files, chromosomes)
       aggregated_files = aggregate(depth_files.groupTuple())
   } else {
       aggregated_files = Channel.fromPath(params.input_files)
   }
   create_accessibility_mask(summarize(aggregated_files))
}