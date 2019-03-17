#' @author Kevin Chau
#' @created 2019-03-09
#' @modified 2019-03-16
#' @description tSNE calculation and visualization of aligned seurat object

library(Seurat)
library(tidyverse)

seurat_align <- readRDS(snakemake@input[[1]])

# TSNE construction

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

print(snakemake@output[[1]])
saveRDS(seurat_tsne, snakemake@output[[1]])

# Visualizations

# Seurat objects for individual genotypes

ident_use <- sort(unique(seurat_tsne@meta.data$name))

seurat_indiv <- setNames(
    lapply(
        ident_use,
        function(i) {
            StashIdent(seurat_tsne, save.name = "clusterID") %>%
                SetAllIdent(., id = "name") %>%
                SubsetData(., ident.use = i) %>%
                SetAllIdent(., id = "clusterID")
        }
    ),
    ident_use
)

tsne_gg <- function(srt, cluster_labels = NULL) {
    pt_data <- cbind(
        srt@dr$tsne@cell.embeddings,
        srt@meta.data[rownames(srt@dr$tsne@cell.embeddings), ]
    ) %>%
        mutate(ident = srt@ident)
    if (!is.null(cluster_labels)) {
        pt_data$label <- cluster_labels[as.character(pt_data$ident)]
    } else {
        pt_data$label <- pt_data$ident
    }
    ret_pt <- ggplot() +
        geom_point(
            data = pt_data, 
            aes(x = tSNE_1, y = tSNE_2, colour = ident), 
            size = 1.5
        ) +
        geom_text(
            data = pt_data %>% 
                group_by(ident) %>% 
                summarise(
                    med_tsne_1 = median(tSNE_1), med_tsne_2 = median(tSNE_2)
                ),
            aes(
                x = med_tsne_1, y = med_tsne_2, label = ident
            ),
            colour = "black", size = 14, fontface = "bold"
        ) +
        guides(colour = guide_legend(title = "")) +
        scale_colour_discrete(
            labels = unique(pull(arrange(pt_data, ident), label))
        )
    return(ret_pt)
}

# Merged plot

merge_pt <- tsne_gg(seurat_tsne, cluster_labels = NULL) +
    ggtitle("MERGE")
ggsave(
    filename = snakemake@output[[2]],
    plot = merge_pt,
    width = 16, height = 9, device = "pdf"
)

# Individual plots
try(dev.off())
pdf(snakemake@output[[3]], width = 16, height = 9)
for (i in ident_use) {
    print(tsne_gg(seurat_indiv[[i]], cluster_labels = NULL))
}
dev.off()
