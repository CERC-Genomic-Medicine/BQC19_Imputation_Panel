# Minimac3 and Minimac4 Imputation Pipeline

This repository provides a Nextflow pipeline for performing imputation using Minimac3 and Minimac4 software. The pipeline enables imputation of genetic data, including chromosome X. First, by generating a composite reference panel using Minimac3 and then by inferring the unobserved genotypes in the genotyping array data using Minimac4.

## Introduction

Imputation plays a crucial role in genotype analysis by inferring unobserved genetic variants based on a reference panel. Minimac3 and Minimac4 are widely used software tools for accurate genotype imputation. This pipeline streamlines the imputation process, allowing for parallel imputation of autosomal chromosomes as well as chromosome X.

## Pipeline Steps

The pipeline follows these main steps:

1. **Data Preparation:** The pipeline requires study samples and a reference panel in VCF format, one vcf per each autosomal chromosomes as well as one VCF for psuedoautosomal region of chromosome X and one VCF for non-psuedoautosomal region of chromosome X for males and one for non-psuedoautosomal region of chromosome X for females(please include term female and male in the name of the vcf files). Don't forget to set chrX parameter of nextflow config file as true if you are planning to run the pipeline for chromosome X. You need to run the pipeline for autosomes and chrX seperately. Make sure to organize your data accordingly before running the pipeline.

2. **Autosomal Imputation:** Perform genotype imputation separately for autosomal chromosomes using Minimac3 or Minimac4 software. This step leverages a reference panel to impute missing variants.

3. **Chromosome X Imputation:** Perform genotype imputation separately for chromosome X using Minimac3 or Minimac4 software. Imputation is performed for males and females seperatly in non-pseudoautosomal regions by accounting for the differences in this region.

4. **Output Generation:** The pipeline generates imputed genotype data in VCF format, as well as empirical dosage files in VCF format, providing a valuable resource for downstream analyses such as association studies or population genetics.

## Input Data

The pipeline expects the following input data:

- Study samples: Genotype data for the study samples in VCF format and their VCF index files per chromosome
- Reference panel: Haplotype reference panel in VCF format and their VCF index files per chromosome

## Output

The pipeline generates the following outputs:

- Imputed genotypes as well as empirical dosage files: The study samples' genotypes are imputed separately for autosomal chromosomes and chromosome X using Minimac4 in VCF format.

## Running the Pipeline

To run the pipeline, follow the steps below:

1. Ensure you have the required dependencies installed, including Minimac3 and Minimac4 and Nextflow.

2. Clone this repository to your local machine or your account on the remote cluster using the following command:

```bash
git clone https://github.com/CERC-Genomic-Medicine/BQC19_Imputation_Panel.git
```
3. Create a new folder where the pipeline results and log files will be generated.

4. Navigate to the new folder and copy the nextflow.config file from this repository.

5. Edit the nextflow.config file to customize it according to your requirements. Provide the file paths for the imputed and truth VCF files. Please note that this pipeline only supports VCF format.

6. Execute the pipeline using the provided command:
```bash
sbatch --account="name of the account" --time=48:00:00 --mem=4G -J imputation --wrap="nextflow run /path/to/imputation.nf" -o imputation.slurm.log
```

## Dependencies

Make sure you have the following dependencies installed:

- [Minimac3](https://genome.sph.umich.edu/wiki/Minimac4)
- [Minimac4](https://genome.sph.umich.edu/wiki/Minimac4)
- [Nextflow](https://www.nextflow.io/)
- [bcftools](https://samtools.github.io/bcftools/howtos/install.html)
- [tabix](https://howtoinstall.co/en/tabix)

Acknowledgements
We would like to acknowledge the developers of Minimac software for providing the tools necessary for imputation. Their contribution is greatly appreciated.