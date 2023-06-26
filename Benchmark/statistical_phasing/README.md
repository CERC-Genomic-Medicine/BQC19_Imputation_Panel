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

1. Ensure you have the required dependencies installed, including Beagle 5.4 software.
2. Prepare the input data in the appropriate format.
3. Execute the pipeline using the following command:

```bash
sbatch --account="name of the account" --time=48:00:00 --mem=4G -J phasing --wrap="nextflow run /path/to/statisticalphasing.nf" -o phasing.slurm.log
```
## Dependencies

Ensure you have the following dependencies installed:

- [Nextflow](https://www.nextflow.io/)
- [BEAGLE software](http://faculty.washington.edu/browning/beagle/beagle.html#download)
- [bcftools](https://samtools.github.io/bcftools/howtos/install.html)

Acknowledgements
We would like to acknowledge the developers of BEAGLE software for providing the tools necessary for phasing and imputation. Their contribution is greatly appreciated.