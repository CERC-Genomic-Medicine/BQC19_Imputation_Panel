process merge_ref {
  cache "lenient"
  memory "32GB"
  time "4h"
  scratch true

  input:
  path (vcfs)

  output:
  path "reference.vcf.gz"

  publishDir "result/merged_ref/", pattern: "*.vcf.gz", mode: "copy"

  script:
  """
  bcftools concat $vcfs -Oz -o reference.vcf.gz
  """
} 

process merge_study {
  cache "lenient"
  memory "32GB"
  time "4h"
  scratch true

  input:
  path (vcfs)

  output:
  path "study.vcf.gz"

  publishDir "result/merged_study/", pattern: "*.vcf.gz", mode: "copy"

  script:
  """
  bcftools concat $vcfs -Oz -o study.vcf.gz
  """
} 

process convert_geno {
  cache "lenient"
  memory "64GB"
  time "2h"
  scratch true

  input:
  path(vcf)
  
  output:
  tuple path("reference.geno"), path("reference.site")

  publishDir "result/reference_geno_site/", pattern: "*.geno", mode: "copy"
  publishDir "result/reference_geno_site/", pattern: "*.site", mode: "copy"
  script:
  """
  ${params.path_to_laser}/vcf2geno/vcf2geno --inVcf $vcf --out reference
  """
}


process reference_PCA {
  cache "lenient"
  memory "32GB"
  time "1h"
  scratch true

  input:
  tuple path(geno), path(site)

  output:
  path("reference.RefPC.coord")
  path(geno)
  path(site)

  publishDir "result/reference_PCA/", pattern: "*.coord", mode: "copy"

  script:
  """
  ${params.path_to_laser}/laser -g $geno -k ${params.nPCs} -pca 1 -o reference
  """

}

process chunk_sample {
    cache "lenient"
    memory "8GB"
    time "1h"
    scratch true


    input:
    path(study_vcf)
    output:
    path("samples.*.txt")
    path(study_vcf)
    
    publishDir "result/samples/", pattern: "*.txt", mode: "copy"

    """
    bcftools  view -h $study_vcf  | tail -n1 | cut -f10- | tr '\t' '\n' > "samples.txt"
    split -l 20 -d --additional-suffix=.txt "samples.txt" samples.
    """
}

process convert_geno2 {
  cache "lenient"
  errorStrategy 'retry'
  maxRetries 3
  memory "8GB"
  time "1h"
  scratch true

  input:
  path(samples)
  each study_vcf

  output:
  tuple path("*.geno"), path("*.site")

  publishDir "result/study_geno_site/", pattern: "*.geno", mode: "copy"
  publishDir "result/study_geno_site/", pattern: "*.site", mode: "copy"

  script:
  """
  bcftools view -S $samples  $study_vcf -Oz -o  ${samples.getBaseName()}.vcf.gz
  ${params.path_to_laser}/vcf2geno/vcf2geno --inVcf ${samples.getBaseName()}.vcf.gz --out ${samples.getBaseName()}
  """

}


process project {
  errorStrategy 'retry'
  maxRetries 3
  memory "32GB"
  time "3h"
  cpus "8"

  input:
  each path(ref_coord)
  each path(reference_geno)
  each path(reference_site)
  tuple path(study_geno), path(study_site)

  output:
  path "*.ProPC.coord"

  publishDir "result/study_PCA/", pattern: "*.coord", mode: "copy"

  script:
  """
  ${params.path_to_laser}/trace -s $study_geno -g $reference_geno -c $ref_coord -k 20 -K 100 -l 1000 -nt 8  -o ${study_geno.getBaseName()}
  """

}


process concat_PC_files {
    cache "lenient"
    cpus 1  
    memory "4GB"
    time "00:15:00"
    scratch true

    input:
    path(proCoord)
    output:     
    path("all_samples_coord.ProPC.coord")
        
    publishDir "result/study_ProPCs_all/", pattern: "*.ProPC.coord", mode: "copy"


   """
   awk 'FNR>1' ${proCoord} > all_samples_coord.ProPC.coord
   sed -i -e '1ipopID\tindivID\tL\tK\tt\tZ\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10\tPC11\tPC12\tPC13\tPC14\tPC15\tPC16\tPC17\tPC18\tPC19\tPC20\n' all_samples_coord.ProPC.coord
   """
   
}

  
process infer_ancestry {
  cache "lenient"
  cpus 1  
  memory "8GB"
  time "00:30:00"
  scratch true

  input:
  path(ref_coord)
  path(study_coord)

  output:
  path "predicted_ancestry.txt"

  publishDir "result/", mode: "copy"

  script:
  """
  ancestry_estimation.py -r $ref_coord -s $study_coord -rl ${params.reference_labels} -np ${params.nPCs}
  """  

}

workflow {
ref_ch = Channel.fromPath(params.input_reference) | map { file -> [ file, file + ".tbi"] }
study_ch = Channel.fromPath(params.input_study) | map { file -> [ file, file + ".tbi"] }

reference_vcf = merge_ref(ref_ch.collect())
study_vcf = merge_study(study_ch.collect())

reference_geno = convert_geno(reference_vcf)
(subset_samples, vcf)  = chunk_sample(study_vcf)

study_geno = convert_geno2(subset_samples.flatten(), vcf)

(reference_PC_coord, ref_geno, ref_site) = reference_PCA(reference_geno)

study_proPC_coord = project(reference_PC_coord, ref_geno, ref_site, study_geno)
study_proPC_coords = concat_PC_files(study_proPC_coord.collect())
infer_ancestry(reference_PC_coord, study_proPC_coords)

}
