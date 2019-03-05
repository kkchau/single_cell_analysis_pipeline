# Single-Cell RNA-seq Analysis with CellRanger and Seurat

Author: Kevin Chau  
Date: 2019-03-03

Snakemake pipeline for single-cell RNA-seq analysis using cellranger and Seurat.

## Usage

Since downstream Seurat analysis requires some custom parameter inputs, this pipeline is divided into several parts:
1. CellRanger is used to quantify the RNA-seq data and an initial Seurat object is initialized per sample.
    * Browse the results of this portion of the pipeline and tune the respective parameters in `config.yaml`; details are outlined in the file
2. This Seurat object is then processed utilizing user-defined parameters in `config.yaml`.

### Initial counts and setup
