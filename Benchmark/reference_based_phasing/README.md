# EAGLE 2 Pre-Phasing Pipeline

This repository contains a Nextflow pipeline for performing pre-phasing of study samples using the EAGLE 2 software. The pipeline utilizes a reference panel and performs pre-phasing separately for autosomal chromosomes as well as chromosome X. Furthermore, it performs phasing separately for males and females in the non-pseudoautosomal parts of chromosome X.

## Introduction

Pre-phasing is an essential step in genotype imputation and haplotype analysis. EAGLE 2 is a widely used software tool for accurate pre-phasing of genetic data. This pipeline streamlines the pre-phasing process by leveraging EAGLE 2 and allowing for parallel pre-phasing of autosomal chromosomes and chromosome X.

## Pipeline Steps

The pipeline follows these main steps:

1. **Data Preparation:** The pipeline requires study samples and a reference panel in VCF format, one vcf per each autosomal chromosomes as well as one VCF for psuedoautosomal region of chromosome X and one VCF for non-psuedoautosomal region of chromosome X for males and one for non-psuedoautosomal region of chromosome X for females(please include term female and male in the name of the vcf files). Don't forget to set chrX parameter of nextflow config file as true if you are planning to run the pipeline for chromosome X. You need to run the pipeline for autosomes and chrX seperately. Make sure to organize your data accordingly before running the pipeline.

2. **Autosomal Pre-Phasing:** The pipeline performs pre-phasing of the study samples separately for autosomal chromosomes using EAGLE 2. This step leverages the reference panel to infer accurate haplotypes.

3. **Chromosome X Pre-Phasing:** The pipeline performs pre-phasing of the study samples separately for chromosome X. It also performs phasing in the non-pseudoautosomal regions of chromosome X seperately for males and females.

4. **Output Generation:** The pipeline generates phased haplotypes for the study samples in VCF format, providing a valuable resource for downstream analyses such as imputation or association studies.

## Input Data

The pipeline expects the following input data:

- Study samples: Genotype data for the study samples in VCF format and their VCF index files per chromosome.
- Reference panel: Haplotype reference panel in VCF format and their VCF index files per chromosome.

## Output

The pipeline generates the following outputs:

- Phased haplotypes: The study samples' haplotypes are phased separately for autosomal chromosomes and chromosome X (pseudoautosomal, and non-pseudoautosomal regions) using EAGLE 2.

## Running the Pipeline

To run the pipeline, follow the steps below:

1. Ensure you have the required dependencies installed, including EAGLE2 and Nextflow.

2. Clone this repository to your local machine or your account on the remote cluster using the following command:

```bash
git clone https://github.com/CERC-Genomic-Medicine/BQC19_Imputation_Panel.git
```
3. Create a new folder where the pipeline results and log files will be generated.

4. Navigate to the new folder and copy the nextflow.config file from this repository.

5. Edit the nextflow.config file to customize it according to your requirements. Provide the file paths for the imputed and truth VCF files. Please note that this pipeline only supports VCF format.

6. Execute the pipeline using the provided command:

```bash
sbatch --account="name of the account" --time=48:00:00 --mem=4G -J prephasing --wrap="nextflow run /path/to/phasing.nf" -o phasing.slurm.log
```
## Dependencies

Make sure you have the following dependencies installed:

- [EAGLE 2](https://alkesgroup.broadinstitute.org/Eagle/) 
- [Nextflow](https://www.nextflow.io/)
- [bcftools](https://samtools.github.io/bcftools/howtos/install.html)
- [tabix](https://howtoinstall.co/en/tabix)


Acknowledgements
We would like to acknowledge the developers of EAGLE software for providing the tools necessary for pre-phasing. Their contribution is greatly appreciated.