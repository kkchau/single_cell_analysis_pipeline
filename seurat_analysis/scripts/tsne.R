#' @author Kevin Chau
#' @created 2019-03-09
#' @modified 2019-03-09
#' @description tSNE calculation and visualization of aligned seurat object

library(Seurat)
library(tidyverse)

seurat_align <- readRDS(snakemake@input[[1]])

seurat_tsne <- RunTSNE(
    object = seurat_align, reduction.use = "cca.aligned",
    dims.use = seq(
        as.numeric(snakemake@params[["align_dims_start"]]),
        as.numeric(snakemake@params[["align_dims_end"]])
    )
) %>%
    FindClusters(
        object = ., reduction.type = "cca.aligned", 
        resolution = as.numeric(snakemake@params[["tsne_resolution"]])
    )

saveRDS(seurat_tsne, snakemake@output[["tsne"]])