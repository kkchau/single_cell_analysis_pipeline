library(Seurat)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

results_dir <- list()
results_dir[args[1]] <- snakemake@input[1]

seurat_data <- Read10X(
    data.dir = results_dir
) %>%
    CreateSeuratObject(min.cells = 3) %>%
    AddMetaData(
        .,
        metadata = (
            Matrix::colSums(.@raw.data[
                grep(pattern = "^mt-|^MT-", rownames(.@raw.data), value = TRUE),
            ]) /
                Matrix::colSums(.@raw.data)
        ),
        col.name = "percent_mito"
    )