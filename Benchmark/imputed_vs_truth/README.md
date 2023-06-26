# Imputed Vs Truth Pipeline

This Nextflow pipeline is designed to compare imputed VCF files of a reference panel with truth VCF files. It performs various processes to analyze and evaluate the concordance between the imputed and truth data.

## Pipeline Processes

The pipeline consists of the following processes:

1. **concat_imputed Process**
This process concat all the imputed VCF files for all chromosomes.

2. **concat_truth Process**
This process concat all the WGS VCF files for all chromosomes.

3. **get_imputed_sample_names Process**
This process extracts the sample names from the imputed VCF files.

4. **imputed_vs_truth Process**
This process runs the `imputed_vs_truth.py` Python script, which iterates over both the truth and imputed VCF files. It assigns a specific label to each variant based on its presence in both truth and imputed files and the number of alternate alleles. The process generates compressed TXT files containing information for each individual, including chromosome number, position, alternate allele, reference allele, imputed genotype, truth genotype, and label.

5. **concat_by_sample Process**
This process concatenates all the TXT files generated in the previous process, resulting in one large TXT file per individual that contains information for all chromosomes.

6. **generate_summary Process**
This process calculates the number of variants in each category and introduces new metrics for concordance calculation. The generated summary includes metrics such as count of variants in whole-genome sequencing (WGS) data, variants present in both imputed and truth data, concordance percentages, coverage metrics, and more.

7. **concat_all_samples_summary Process**
This process concatenates all the summary files for all individuals, creating a single TXT file.

## Input Data

The pipeline expects the following input data to evaluate the performance of reference panel for imputation:

- Imputation files for the evaluation set: Imputed files in VCF format as well as VCF index files. 
- Whole genome sequencing files for the evaluation set:  Whole genome sequencing files in VCF format as well as VCF index files. 

## Output

The pipeline generates the following outputs:

- Imputation quality files per individual in the evaluation set: These files contains list of variants present in the imputed data and WGS data and a label for each variant showing whether the variant is concordantly imputed or not.
- Summary file: Contains metrics such as count of variants in whole-genome sequencing (WGS) data, variants present in both imputed and truth data, concordance percentages, coverage metrics, and more.


## Running the Pipeline

1. Ensure you have the required dependencies installed, including Nextflow.

2. Clone this repository to your local machine or your account on the remote cluster using the following command:

```bash
git clone https://github.com/CERC-Genomic-Medicine/BQC19_Imputation_Panel.git
```
3. Create a new folder where the pipeline results and log files will be generated.

4. Navigate to the new folder and copy the nextflow.config file from this repository.

5. Edit the nextflow.config file to customize it according to your requirements. Provide the file paths for the imputed and truth VCF files. Please note that this pipeline only supports VCF format.

6. Modify and execute the following bash commands in your folder. Alternatively, you can generate a bash file using these commands and then execute it:

```bash 
#!/bin/bash

module load tabix
module load nextflow

source /path/to/your/python/virtual/env/bin/activate
sbatch --account=account_name --time=24:00:00 --mem=16G -J Post_Imputation --wrap="nextflow run /path/to/BQC19_Imputation_Panel/Benchmark/imputed_vs_truth/ImputedVsTruth.nf" -o post_imp.slurm.log
deactivate
```

## Dependencies

- Python packages: pysam, pandas
- [Nextflow](https://www.nextflow.io/)
- [bcftools](https://samtools.github.io/bcftools/howtos/install.html)
- [tabix](https://howtoinstall.co/en/tabix)


