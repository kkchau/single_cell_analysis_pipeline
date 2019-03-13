#' @author Kevin Chau
#' @created 2019-03-05
#' @modified 2019-03-08
#' @description Align the subspaces of the merged Seurat object

library(Seurat)
library(tidyverse)

seurat_merge <- readRDS(snakemake@input[[1]])

# Visualze CCA

pdf(snakemake@output[[1]], width = 16, height = 9)

DimPlot(
    object = seurat_merge, 
    reduction.use = "cca", group.by = "name", 
    pt.size = 0.5, do.return = TRUE
)
VlnPlot(
    object = seurat_merge, 
    features.plot = "CC1", group.by = "name", 
    do.return = TRUE
)

MetageneBicorPlot(
    seurat_merge, grouping.var = "name", 
    dims.eval = 1:30, display.progress = TRUE
)

DimHeatmap(
    object = seurat_merge, reduction.type = "cca", 
    cells.use = 500, 
    dim.use = seq(
        snakemake@params[["heatmap_start"]],
        snakemake@params[["heatmap_end"]]
    ), 
    do.balanced = TRUE, 
    do.return = TRUE
)

# Align subspaces

seurat_align <- AlignSubspace(
    seurat_merge, 
    reduction.type = "cca", 
    grouping.var = "name", 
    dims.align = seq(
        snakemake@params[["align_dims_start"]],
        snakemake@params[["align_dims_end"]]
    )
)
saveRDS(seurat_align, snakemake@output[["cca_align"]])