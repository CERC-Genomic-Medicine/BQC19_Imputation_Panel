#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>
* Adapted and modified get_chunk and ligate from https://github.com/CERC-Genomic-Medicine/shapeit4_pipeline/blob/main/Phasing.nf by Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 3.0
* YEAR: 2023
*/

process get_ref_chr_names{
    cache "lenient"
    cpus 1
    memory "4GB"
    time "00:30:00"
    scratch true

    input:
    tuple val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple stdout, val(sex_id), path(vcf), path(vcf_index)

    """
    tabix -l ${vcf} 
    """
}

process rm_chr_name_ref{
    cache "lenient"
    cpus 1
    memory "4GB"
    time "00:30:00"

    scratch true
    input:
    tuple val(chr_name), val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple val(chr_name), val(sex_id), path("*.vcf.gz"), path("*.vcf.gz.tbi")
    
    script:
    if (params.chrX == false) {
    """
    echo "${chr_name}" > chroms1.txt
    chr=${chr_name}
    echo "\${chr:3}" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.vcf.gz
    """
    } else {
    """
    echo "${chr_name}" > chroms1.txt
    echo "X" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.vcf.gz
    """
    }
}


process get_study_chr_names{
    cache "lenient"
    cpus 1
    memory "4GB"
    time "00:30:00"
    scratch true

    input:
    tuple val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple stdout,  val(sex_id), path(vcf), path(vcf_index)

    """
    tabix -l ${vcf}
    """
}

process rm_chr_name_study{
    cache "lenient"
    cpus 1
    memory "4GB"
    time "00:30:00"
    scratch true  

    input:
    tuple val(chr_name), val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple val(chr_name), val(sex_id), path("*.vcf.gz"), path("*.vcf.gz.tbi")
    
    script:
    if (params.chrX == false) {
    """
    echo "${chr_name}" > chroms1.txt
    chr=${chr_name}
    echo "\${chr:3}" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.vcf.gz
    """
    } else {
    """
    echo "${chr_name}" > chroms1.txt
    echo "X" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.vcf.gz
    """
    }
}


process convert_ref_vcf{
    //errorStrategy 'retry'
    //maxRetries 3
    cache "lenient"
    cpus 4
    memory "64GB"
    time "5h"
    scratch true

    input:
    tuple val(chr_name), val(sex_id), path(ref_vcf), path(ref_vcf_index)
    
    output:
    tuple val(chr_name), val(sex_id), path( "*.m3vcf.gz")
    
    publishDir "minimac_m3vcfs/", pattern: "*.m3vcf.gz", mode: "copy"

    """
    ${params.minimac3} --refHaps ${ref_vcf} --cpus 4 --processReference --prefix ${ref_vcf.getBaseName()}
    """
}

process get_chunks_study {
	//executor "local"
	cache "lenient"
	cpus 1
	memory "4GB"
    time "00:30:00"

	input:
	tuple val(chr_name), val(start_bp), val(stop_bp), val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple path("*.chunk"), val(chr_name), val(sex_id), path(vcf), path(vcf_index)

	"""
    chrom=`bcftools index -s ${vcf} | cut -f1`
       
    extend=0
	for i in `seq ${start_bp} ${params.window} ${stop_bp}`; do
		if [ \${extend} -eq 0 ]; then
			chunk_start=\$((${params.flank} > i ? 1 : i - ${params.flank}))
		fi
		chunk_stop=\$((i + ${params.window} + ${params.flank}))
		n=`bcftools view -HG ${vcf} \${chrom}:\${chunk_start}-\${chunk_stop} | wc -l`
		if [ \${n} -gt 0 ]; then
			printf "\${chrom}\t\${chunk_start}\t\${chunk_stop}\t\${n}\n" > \${chrom}_\${chunk_start}_\${chunk_stop}.chunk
			extend=0
		else
			extend=1
		fi
	done
	if [ \${extend} -eq 1 ]; then
		printf "\${chrom}\t\${chunk_start}\t\${chunk_stop}\t\${n}\n" > \${chrom}_\${chunk_start}_\${chunk_stop}.chunk
	fi
	"""
}

process minimac_imputation{
    //errorStrategy 'retry'
    //maxRetries 3
    cache "lenient"
    cpus 4
    memory "128GB"
    time "10h"
    scratch true  

    input:
    tuple val(chr_name), val(sex_id), file(ref_vcf), path(chunk), file(study_vcf), file(study_vcf_index)
     
    output:
    tuple val(chr_name), path("*.imp.dose.vcf.gz"), path("*.imp.dose.vcf.gz.tbi"), path("*.imp.empiricalDose.vcf.gz"), path("*.imp.empiricalDose.vcf.gz.tbi"), path("*.info")

    publishDir "imputed_dose_vcfs/", pattern: "*.imp.dose.vcf.gz*", mode: "copy"
    publishDir "empirical_dose_vcfs/", pattern: "*.imp.empiricalDose.vcf.gz*", mode: "copy"
    publishDir "imputed_info/", pattern: "*.info", mode: "copy"

    script:
    if (params.chrX == false) {
    """
    chrom=`head -n1 ${chunk} | cut -f1`
    start_bp=`head -n1 ${chunk} | cut -f2`
	stop_bp=`head -n1 ${chunk} | cut -f3`
    bcftools view -r \${chrom}:\${start_bp}-\${stop_bp} ${study_vcf} -Oz -o study.\${chrom}_\${start_bp}_\${stop_bp}.vcf.gz 
    bcftools index --tbi study.\${chrom}_\${start_bp}_\${stop_bp}.vcf.gz

    ${params.minimac4} --refHaps $ref_vcf --haps study.\${chrom}_\${start_bp}_\${stop_bp}.vcf.gz  --chr \${chrom} --start \${start_bp} --end \${stop_bp} --minRatio 0.00001 --prefix study.\${chrom}_\${start_bp}_\${stop_bp} --cpus 4  --meta --ignoreDuplicates

    echo "${chr_name}" > chroms1.txt
    chr=${chr_name}
    echo "\${chr:3}" > chroms2.txt

    paste chroms2.txt chroms1.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt study.\${chrom}_\${start_bp}_\${stop_bp}.dose.vcf.gz -Oz -o study.\${chrom}_\${start_bp}_\${stop_bp}.imp.dose.vcf.gz
    bcftools index --tbi study.\${chrom}_\${start_bp}_\${stop_bp}.imp.dose.vcf.gz

    bcftools annotate --rename-chrs chr_name_conv.txt study.\${chrom}_\${start_bp}_\${stop_bp}.empiricalDose.vcf.gz -Oz -o study.\${chrom}_\${start_bp}_\${stop_bp}.imp.empiricalDose.vcf.gz
    bcftools index --tbi study.\${chrom}_\${start_bp}_\${stop_bp}.imp.empiricalDose.vcf.gz
    """
    } else {
    """
    chrom=`head -n1 ${chunk} | cut -f1`
    start_bp=`head -n1 ${chunk} | cut -f2`
	stop_bp=`head -n1 ${chunk} | cut -f3`
    bcftools view -r \${chrom}:\${start_bp}-\${stop_bp} ${study_vcf} -Oz -o study.\${chrom}_\${start_bp}_\${stop_bp}.vcf.gz 
    bcftools index --tbi study.\${chrom}_\${start_bp}_\${stop_bp}.vcf.gz

    ${params.minimac4} --refHaps $ref_vcf --haps study.\${chrom}_\${start_bp}_\${stop_bp}.vcf.gz  --chr \${chrom} --start \${start_bp} --end \${stop_bp} --minRatio 0.00001 --prefix study.\${chrom}_\${start_bp}_\${stop_bp} --cpus 4  --meta --ignoreDuplicates

    echo "chrX" > chroms1.txt
    echo "X" > chroms2.txt

    paste chroms2.txt chroms1.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt study.\${chrom}_\${start_bp}_\${stop_bp}.dose.vcf.gz -Oz -o study.\${chrom}_\${start_bp}_\${stop_bp}.imp.dose.vcf.gz
    bcftools index --tbi study.\${chrom}_\${start_bp}_\${stop_bp}.imp.dose.vcf.gz

    bcftools annotate --rename-chrs chr_name_conv.txt study.\${chrom}_\${start_bp}_\${stop_bp}.empiricalDose.vcf.gz -Oz -o study.\${chrom}_\${start_bp}_\${stop_bp}.imp.empiricalDose.vcf.gz
    bcftools index --tbi study.\${chrom}_\${start_bp}_\${stop_bp}.imp.empiricalDose.vcf.gz
    """
    }
}

process ligate {
	
	cache "lenient"
	scratch true
	cpus 1
	memory "32G"
	time "12h"

	input:
	tuple val(chromosome), path(imputed_vcfs), path(imputed_vcfs_index), path(imputed_emp_vcfs), path(imputed_emp_vcfs_index), path(info)

	output:
	tuple path("*.imputed.dose.vcf.gz"), path("*.imputed.dose.vcf.gz.tbi"), path("*.imputed.empiricalDose.vcf.gz"), path("*.imputed.empiricalDose.vcf.gz.tbi")

	publishDir "final_imputed_vcfs/", pattern: "*.vcf.gz*", mode: "copy"

	"""
	for f in ${imputed_vcfs}; do echo \${f}; done | sort -V > files_list.txt
	bcftools concat -f files_list.txt -l -Oz -o ${chromosome}.imputed.dose.vcf.gz
	bcftools index --tbi ${chromosome}.imputed.dose.vcf.gz
	for f in ${imputed_emp_vcfs}; do echo \${f}; done | sort -V > files_list.txt
    bcftools concat -f files_list.txt -l -Oz -o ${chromosome}.imputed.empiricalDose.vcf.gz
    bcftools index --tbi ${chromosome}.imputed.empiricalDose.vcf.gz
	"""
}


workflow {

        ref_ch = Channel.fromPath(params.ref_vcf_path).map{ vcf -> [ vcf.name.toString().tokenize('.').contains('female'), vcf, vcf + ".tbi" ] }
        study_ch = Channel.fromPath(params.study_vcf_path).map{ vcf -> [ vcf.name.toString().tokenize('.').contains('female'), vcf, vcf + ".tbi" ] }

        chromosome_sizes = Channel.fromList([['chr1\n', 10352, 248946390], ['chr2\n', 10369, 242183307], ['chr3\n', 10374, 198230061], ['chr4\n', 10035, 190204491], ['chr5\n', 11507, 181477936], ['chr6\n', 63746, 170744468], ['chr7\n', 16359, 159335877],
    ['chr8\n', 72806, 145075995], ['chr9\n', 10421, 138318560], ['chr10\n', 10331, 133787400], ['chr11\n', 70402, 135076530], ['chr12\n', 10241, 133264323], ['chr13\n', 16000260, 114353949], ['chr14\n', 16022732, 106883658], ['chr15\n', 17000328, 101981149],
    ['chr16\n', 10085, 90228306], ['chr17\n', 104112, 83244562], ['chr18\n', 10652, 80262913], ['chr19\n', 60603, 58605828], ['chr20\n', 60148, 64333628], ['chr21\n', 5030957, 46699945], ['chr22\n', 10519389, 50808269]])

        ref_vcfs = get_ref_chr_names(ref_ch)
        study_vcfs = get_study_chr_names(study_ch)
        
        ref_rm_chr_vcfs = rm_chr_name_ref(ref_vcfs)
        study_rm_chr_vcfs = rm_chr_name_study(study_vcfs)

        study_chunk_ch = chromosome_sizes.join(study_rm_chr_vcfs, by:[0])
        //ref_chunks = get_chunks_ref(ref_rm_chr_vcfs)
        //chunks_all_ref = ref_chunks.flatMap { chunks, chromosome, sex_id, vcf, vcf_index ->
        //chunks.collect { chunk -> [chunk, chromosome, sex_id, vcf, vcf_index] }
        //}
        ref_cnv_vcfs = convert_ref_vcf(ref_rm_chr_vcfs)

        study_chunks = get_chunks_study(study_chunk_ch)
        chunks_all_study = study_chunks.flatMap { chunks, chromosome, sex_id, vcf, vcf_index ->
        chunks.collect { chunk -> [chromosome, sex_id, chunk, vcf, vcf_index] }
        }
        imputation_ch = ref_cnv_vcfs.combine(chunks_all_study, by:[0, 1])
        imputed_chunks = minimac_imputation(imputation_ch)
        ligate(imputed_chunks.groupTuple())

}