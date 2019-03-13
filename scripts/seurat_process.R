#' @author Kevin Chau
#' @created 2019-03-04
#' @modified 2019-03-04
#' @description Process Seurat object

library(Seurat)
library(tidyverse)

print(snakemake@input[[1]])

name <- snakemake@wildcards[["sample_id"]]
seurat_data <- readRDS(snakemake@input[[1]]) %>%
    SetAllIdent(., id = "name")

low_thresholds <- c(
    as.numeric(snakemake@params[["nUMI_low"]]), 
    as.numeric(snakemake@params[["percent_mito_low"]])
)
high_thresholds <- c(
    as.numeric(snakemake@params[["nUMI_high"]]), 
    as.numeric(snakemake@params[["percent_mito_high"]])
)

seurat_filter <- FilterCells(
    seurat_data,
    subset.names = c("nUMI", "percent_mito"),
    low.thresholds = c(
        as.numeric(snakemake@params[["nUMI_low"]]), 
        as.numeric(snakemake@params[["percent_mito_low"]])
    ),
    high.thresholds = c(
        as.numeric(snakemake@params[["nUMI_high"]]),
        as.numeric(snakemake@params[["percent_mito_high"]])
    ),
    cells.use = WhichCells(seurat_data, ident = name)
)

# Normalize data, scale, and find variable genes

seurat_scale <- NormalizeData(seurat_filter) %>%
    FindVariableGenes(., do.plot = FALSE) %>%
    ScaleData(., vars.to.regress = c("nUMI", "percent_mito"))


saveRDS(seurat_scale, snakemake@output[[1]])
