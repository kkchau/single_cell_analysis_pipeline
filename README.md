# Single-Cell RNA-seq Analysis with CellRanger and Seurat

Author: Kevin Chau
Date: 2019-03-03

Snakemake pipeline for single-cell RNA-seq analysis using cellranger and Seurat.

## Instructions

### Cellranger Counts (rule counts)

1. Modify `config.yaml` as necessary for samples (see example strucuture under `config['cellranger']['ids']`
    * If running on a cluster, modify cluster rules as necessary in `cluster.yaml`
    * Please see cellranger [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input) for full description of parameters. This pipeline is based on FASTQ structure based on "Scenario: My FASTQs are in an output folder from mkfastq or bcl2fastq, but there are multiple folders per sample index, like 'SI-GA-A1_1' and 'SI-GA-A1_2'"

### Seurat object construction