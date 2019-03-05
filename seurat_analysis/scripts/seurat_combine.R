#' @author Kevin Chau
#' @created 2019-03-04
#' @modified 2019-03-04
#' @description Initial merging and processing all Seurat objects

library(Seurat)
library(tidyverse)

seurat_objects <- lapply(snakemake@input, readRDS)

# Gene selection for CCA

hvg <- unlist(sapply(
    seurat_objects, 
    function(so) so@var.genes
)) %>%
    unique() %>%
    intersect(
        .,
        sapply(seurat_objects, function(so) rownames(so@scale.data)) %>%
            unlist() %>%
            unique()
    )

# CCA

if (length(seurat_objects) == 2) {
    integrated <- RunCCA(
        seurat_objects[[1]], seurat_objects[[2]], genes.use = hvg, num.cc = 30
    )
} else {
    integrated <- RunMultiCCA(
       seurat_objects, genes.use = hvg, num.cc = 30
    )
}
saveRDS(integrated, snakemake@output[[1]])