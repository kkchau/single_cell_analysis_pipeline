# Single-Cell RNA-seq Analysis with CellRanger and Seurat

Author: Kevin Chau

Snakemake pipeline for single-cell RNA-seq analysis using cellranger and Seurat.

## Usage

Included is a convenient `conda` environment description for ease-of-use and replicability. Please call

```sh
conda env create -f conda.yaml
source activate seurat
```

to initialize this environment prior to running the pipeline

Since downstream Seurat analysis requires some custom parameter inputs, this pipeline is divided into several parts:

1. Subworkflow: **data_setup** 
    * CellRanger is used to quantify the RNA-seq data and an initial Seurat object is initialized per sample.
    * Browse the results of this portion of the pipeline and tune the respective parameters in `data_setup/config.yaml`; details are outlined in the file
2. **Main Workflow**
    * This Seurat object is then processed utilizing user-defined parameters in `config.yaml`.
    * Many scripts from `scripts` are replicated in `user_analysis`. This directory is provided as a convenient workspace for the user.

### Initial counts and setup (**data_setup**)

Run this pipeline first to generate initial Seurat object. User should inspect `data/{sample_id}/violin_and_scatter.pdf` and set `nUMI` and `percent_mito` in `config.yaml` as necessary.

### Automated Seurat analysis (**main**)

After modifying `config.yaml`, the user may run this pipeline to integrate multiple single-cell analyses and perform preliminary analyses.