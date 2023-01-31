# Imputed Vs Truth Pipeline

This nextflow pipeline aims to compare imputed vcf files of a reference panel with truth vcf files. 

The pipeline contains different processes:

### get_imputed_chr_names process
This process extracts the chromosome number from the imputed vcf files.
### get_truth_chr_names process
This process extracts the chromosome number from the truth vcf files.
### get_imputed_sample_names process
This process extracts the sample names from imputed vcf files.
### imputed_vs_truth process
This process runs imputed_vs_truth.py python file, which iterate over both truth and imputed vcf files and aligns a specific label to each variant based on whether that variant is present in both truth and imputed vcf files or not, and if it was present in both, what would be the number of alternate alleles (whether this number is equal in both truth and imputed vcf files.)
Here is the list of labels and their descriptions:
- REF :  variant is only present in imputed files, non-zero genotype [(1, 0), (1, 1), (0, 1)].
- REF_0ALT : variant is only present in imputed files, zero genotype (0, 0).
- WGS_AND_REF_EQ : vartiant is present in both imputed and truth files, non-zero genotype [(1, 0), (1, 1), (0, 1)], number of ALT alleles is equal in imputed vs truth.
- WGS_0ALT_AND_REF_EQ : variant is present in both imputed and truth files, zero genotype (0, 0), number of ALT alleles is equal in imputed vs truth which is equal to zero in this case.
- WGS_AND_REF_LT : variant is present in both imputed and truth files, non-zero genotype [(1, 0), (1, 1), (0, 1)] for truth, number of ALT alleles in truth is less than number of ALT alleles in imputed files.
- WGS_0ALT_AND_REF_LT : variant is present in both imputed and truth files, non-zero genotype [(0, 0)] for truth, number of ALT alleles in truth is less than number of ALT alleles in imputed files (this means that number of alth allels in imputed files for this variant is non-zero)
- WGS_AND_REF_GT : variant is present in both imputed and truth files, non-zero genotype [(1, 0), (1, 1), (0, 1)] for truth, number of ALT alleles in truth is greater than number of ALT alleles in imputed files.
- WGS : variant is only present in truth files, non-zero genotype [(1, 0), (1, 1), (0, 1)].
- WGS_0ALT : variant is only present in truth files, zero genotype [(0, 0)].
The output of the process is compressed txt files containing : ["CHROM", "POS", "ALT", "REF", "IMP_gt", "TRUTH_gt", "label"] for each individual per chromosome.
### concat_by_sample process
This process concats all the txt files from previous process and generates one large txt file per individuals containing information for all chromosomes.
### generate_summary process
This process calculate the number of variants in each category as well as introduced new metrics for concordance calculation. 
Here is the list of categories and metrics:
- Sample_ID : ID of the sample
- WGS : count variants which is 0/1, 1/0, 1/1 in that sample in WGS data.
- WGS_AND_REF : variant in "WGS" and in reference panel.
- WGS_AND_REF_EQ : variant in "WGS_AND_REF" and GT_WGS == GT_REF.
- WGS_AND_REF_LT : variant in "WGS_AND_REF" and GT_WGS < GT_REF.
- WGS_AND_REF_GT : variant in "WGS_AND_REF" and GT_WGS > GT_REF.
- REF : variant in "WGS_AND_REF" and GT_WGS > GT_REF.
- REF_0ALT : variant is in REF and GT_REF == 0/0.
- WGS_0ALT : GT_WGS == 0/0.
- AA_Concordance_PERC : WGS_AND_REF_EQ / WGS_AND_REF.
- AA_Concordance : WGS_AND_REF_EQ
- Coverage : abs(WGS - WGS_AND_REF) / WGS
- RA_Discordance_PERC : (REF - REF_0ALT) / REF
- (REF - REF_0ALT) / REF : (REF - REF_0ALT
### concat_all_samples_summary process
This process concats all summary files for all individuals and returns one txt file

## Usage
1- Create a python virtual environment.
2- Activate your virtual environmnet.
3- Intall following packages in the virtual environment:
- pysam
- pandas
4- Install nextflow.
5- Install tabix.
6- Load nextflow.
7- Load tabix.
8- Make a copy of this repository in your local or in your account on the remote cluster using the following command:
```bash
git clone https://github.com/CERC-Genomic-Medicine/BQC19_Imputation_Panel.git
```
9- Create a new folder in which the results of the pipeline and log files related to running process will be generated.
10- Go to your new folder and copy the nextflow.config file from this repository, there.
11- Edit the nextflow.config file based on your need. Put the address of the imputed and truth vcf files in the nextflow.config file. This pipeline only works with vcf format.
12- Now modify and run the following bash commands in your folder. You can also generate a bash file by following commands and execute that:
```bash
#!/bin/bash

module load tabix
module load nextflow

source /path/to/your/python/virtual/env/bin/activate
sbatch --account=account_name --time=24:00:00 --mem=16G -J Post_Imputation --wrap="nextflow run /path/to/BQC19_Imputation_Panel/Benchmark/imputed_vs_truth/ImputedVsTruth.nf" -o post_imp.slurm.log
deactivate
``` 