# DNA Filtration and Kinship Pipeline

This repository contains a Nextflow DSL2 pipeline for DNA filtration and kinship analysis. The pipeline performs various steps including reading a list of chromosomes, extracting and indexing VCF files, calculating linkage disequilibrium, thinning VCF files, and creating GDS objects for SNP PCA analysis.

## Pipeline Overview
The pipeline performs the following steps:
1. Reads a list of chromosomes.
2. Concatenates all non-thinned VCFs.
3. Runs King on extracted full VCF.
4. For each chromosome, extracts and indexes a VCF using bcftools.
5. Calculates linkage disequilibrium and prunes the VCF using plink2.
6. Thins VCF using plink2 --bp-space 250.
7. Concatenates all thinned VCFs into one final VCF.
8. Creates a GDS object from the final thinned VCF.
9. Creates SNP PCA and outputs covariates.

## Installation
To run this pipeline, you need to have Nextflow installed. You can install Nextflow using the following command:
You also need to have the following dependencies installed:
- bcftools/1.17
- plink/2.00a20230303
- King/2.3.2
- R/4.3.1
- Python/3.11.4
- htslib/1.17

## Usage
To run the pipeline, use the following command:
```sh
nextflow run main_DNA_mega_new.nf -profile slurm
```
## Main Scripts
- `main_DNA_mega_new.nf`: The main Nextflow script that defines the workflow.
- `modules/`: Directory containing Nextflow modules for each step of the pipeline.

## Configuration
The pipeline can be configured using a `nextflow.config` file. You can specify parameters such as input files, output directories, and resource requirements.

## Contributing
Contributions are welcome! Please open an issue or submit a pull request on GitHub.

## License
This project is licensed under the MIT License.
