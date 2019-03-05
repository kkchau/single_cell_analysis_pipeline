#' @author Kevin Chau
#' @created 2019-03-03
#' @modified 2019-03-04
#' @description Initialize Seurat object


library(Seurat)
library(tidyverse)

print(snakemake@wildcards[["sample_id"]])
print(snakemake@params[["data_dir"]])

seurat_data <- Read10X(
    data.dir = as.character(snakemake@params[["data_dir"]])
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
    ) %>%
    AddMetaData(
        ., 
        metadata = data.frame(
            name = snakemake@wildcards[["sample_id"]],
            cell = rownames(.@meta.data),
            row.names = 2
        ), 
        col.name = "name"
    ) %>%
    RenameCells(., add.cell.id = snakemake@wildcards[["sample_id"]])

saveRDS(seurat_data, file.path("data", snakemake@wildcards[["sample_id"]], "seurat.rds"))

# Plotting

dir.create(file.path("data", snakemake@wildcards[["sample_id"]], "figures"))

pdf(file.path("data", snakemake@wildcards[["sample_id"]], "figures", "violin_and_scatter.pdf"), width = 16, height = 9)

VlnPlot(
    seurat_data,
    features.plot = c("nGene", "nUMI", "percent_mito"),
    nCol = 3
)

GenePlot(
    seurat_data,
    gene1 = "nUMI",
    gene2 = "nGene"
)
GenePlot(
    seurat_data,
    gene1 = "nUMI",
    gene2 = "percent_mito",
)

dev.off()
