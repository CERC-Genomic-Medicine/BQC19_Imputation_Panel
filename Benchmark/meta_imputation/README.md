# MetaImputation Pipeline using MetaMinimac

This repository provides a Nextflow pipeline for performing meta-imputation using MetaMinimac software. The pipeline enables the imputation of genetic data by combining imputation results from multiple reference panels, improving the accuracy and coverage of genotype data for downstream analyses.

## Introduction

MetaImputation is an approach that combines imputation results from multiple reference panels to enhance the imputation accuracy and impute rare variants with higher confidence. MetaMinimac is a widely used software tool that implements the MetaImputation framework. This pipeline streamlines the meta-imputation process, allowing for efficient integration of imputation results from two different reference panels.

## Pipeline Steps

The pipeline follows these main steps:

1. **Data Preparation:** Ensure imputation files from both reference panels are prepared in VCF formats, one VCF per chromosome. This pipeline needs dosage and empirical dosage files in the following naming format REFERENCE_NAME.CHROMOSOME_NAME.dose.vcf.gz as well as empirical dosage files in the following naming format REFERENCE_NAME.CHROMOSOME_NAME.empiricalDose.vcf.gz 

2. **Meta-imputation using MetaMinimac:**  This step combines the imputation results from multiple reference panels to generate a more accurate imputed dataset.

3. **Meta-Imputation:** The pipeline integrates the imputed data from each reference panel using MetaMinimac's meta-imputation algorithm. This step leverages the strengths of each reference panel to improve the imputation accuracy, especially for rare variants.

4. **Output Generation:** The pipeline generates the meta-imputed genotype data in VCF format, providing a comprehensive resource for downstream analyses, such as association studies or population genetics.

## Input Data

The pipeline expects the following input data:

- Imputation files from the first reference panel: Imputed dosages and empirical dosage files data using the first panel in VCF format. 
- Imputation files from the second reference panel: Imputed dosages and empirical dosage files data using the second panel  in VCF format. 

Please ensure that VCF file names are in the format that was described earlier in this documentation.

## Output

The pipeline generates the following outputs:

- Meta-imputed genotypes: The study samples' genotypes are meta-imputed using MetaMinimac, combining imputation results from multiple reference panels.

## Running the Pipeline

1. Ensure you have the required dependencies installed, including Metaminimac 2 and Nextflow.

2. Clone this repository to your local machine or your account on the remote cluster using the following command:

```bash
git clone https://github.com/CERC-Genomic-Medicine/BQC19_Imputation_Panel.git
```
3. Create a new folder where the pipeline results and log files will be generated.

4. Navigate to the new folder and copy the nextflow.config file from this repository.

5. Edit the nextflow.config file to customize it according to your requirements. Provide the file paths for the imputed and truth VCF files. Please note that this pipeline only supports VCF format.

6. Execute the pipeline using the provided command:
```bash
sbatch --account="name of the account" --time=48:00:00 --mem=4G -J metaimputation --wrap="nextflow run /path/to/metaimputation.nf" -o metaimputation.slurm.log
```
## Dependencies

Make sure you have the following dependencies installed:

- [MetaMinimac](https://github.com/yukt/MetaMinimac2/tree/master)
- [Nextflow](https://www.nextflow.io/)

Acknowledgements
We would like to acknowledge the developers of Metaminimac software for providing the tools necessary for metaimputation. Their contribution is greatly appreciated.