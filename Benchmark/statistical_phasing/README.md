# Haplotype Reference Panel Construction Pipeline

This repository contains a Nextflow pipeline for constructing a haplotype reference panel for genotype imputation. The pipeline performs statistical haplotype phasing in 25 cM windows for autosomal and X chromosomes independently using Beagle 5.4 software. Additionally, it removes singletons to improve the quality of the reference panel.

## Pipeline Steps

1. **Haplotype Phasing:** In this step, the pipeline performs statistical haplotype phasing separately for autosomal and X chromosomes using Beagle 5.4 software. The phasing is done in 25 cM windows to ensure accurate haplotype construction.

2. **Singleton Removal:** After haplotype phasing, the pipeline removes singletons from the reference panel. Singletons are variants that have only one occurrence in the dataset. Removing singletons helps improve the quality of the reference panel by reducing the chance of including rare or erroneous variants.

## Input Data

The pipeline requires input data in vcf format per chromosome for reference set, as well as vcf index files. 

## Output

The pipeline generates the following outputs:

- **Haplotype Reference Panel:** The constructed haplotype reference panel, consisting of phased haplotypes for autosomal and X chromosomes.

## Running the Pipeline

To run the pipeline, please follow these instructions:

1. Ensure you have the required dependencies installed, including BEAGLE 5.4 and Nextflow.

2. Clone this repository to your local machine or your account on the remote cluster using the following command:

```bash
git clone https://github.com/CERC-Genomic-Medicine/BQC19_Imputation_Panel.git
```
3. Create a new folder where the pipeline results and log files will be generated.

4. Navigate to the new folder and copy the nextflow.config file from this repository.

5. Edit the nextflow.config file to customize it according to your requirements. Provide the file paths for the imputed and truth VCF files. Please note that this pipeline only supports VCF format.

6. Execute the pipeline using the provided command:

```bash
sbatch --account="name of the account" --time=48:00:00 --mem=4G -J phasing --wrap="nextflow run /path/to/statisticalphasing.nf" -o phasing.slurm.log
```
## Dependencies

Ensure you have the following dependencies installed:

- [Nextflow](https://www.nextflow.io/)
- [BEAGLE software](http://faculty.washington.edu/browning/beagle/beagle.html#download)
- [bcftools](https://samtools.github.io/bcftools/howtos/install.html)

Acknowledgements
We would like to acknowledge the developers of BEAGLE software for providing the tools necessary for reference-free phasing. Their contribution is greatly appreciated.