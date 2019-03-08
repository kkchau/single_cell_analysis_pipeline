# Single-Cell RNA-seq Analysis with CellRanger and Seurat

Author: Kevin Chau
Date: 2019-03-04

Snakemake pipeline for single-cell RNA-seq analysis using cellranger and Seurat.

## Usage

Included is a convenient `conda` environment description for ease-of-use and replicability. Please call

```sh
conda create -f conda.yaml
source activate seurat
```

to initialize this environment prior to running the pipeline

NOTE: SLURM submission scripts include the `source` environment call before running snakemake.

Since downstream Seurat analysis requires some custom parameter inputs, this pipeline is divided into several parts:

1. CellRanger is used to quantify the RNA-seq data and an initial Seurat object is initialized per sample.
    * Browse the results of this portion of the pipeline and tune the respective parameters in `data/data_setup/config.yaml`; details are outlined in the file
2. This Seurat object is then processed utilizing user-defined parameters in `seurat_analysis/config.yaml`.
3. Many scripts from `seurat_analysis/scripts` are replicated in `user_analysis`. This directory is provided as a convenient workspace for the user.

### Initial counts and setup (data_setup)

Run this pipeline first to generate initial Seurat object. User should inspect `data/{sample_id}/violin_and_scatter.pdf` and set `nUMI` and `percent_mito` in `seurat_analysis/config.yaml` as necessary.

### Automated Seurat analysis (seurat_analysis)

After modifying `seurat_analysis/config.yaml`, the user may run this pipeline to integrate multiple single-cell analyses and perform preliminary analyses.