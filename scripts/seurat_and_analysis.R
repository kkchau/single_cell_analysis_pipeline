# XH1: HET
# XH2: WT

dir.create("../data/figures")

library(Seurat)
library(tidyverse)
library(ggsignif)

seurat_data <- Read10X(
    data.dir = c(
        HET = "../data/counts_results/XH1/outs/filtered_feature_bc_matrix/",
        WT = "../data/counts_results/XH2/outs/filtered_feature_bc_matrix/"
    )
) %>%
    CreateSeuratObject(project = "Cul3Mice", min.cells = 3) %>%
    AddMetaData(
        .,
        metadata = (
            Matrix::colSums(.@raw.data[
                grep(pattern = "^mt-", rownames(.@raw.data), value = TRUE), 
            ]) /
                Matrix::colSums(.@raw.data)
        ), 
        col.name = "percent_mito"
    )

# Plots

umi_vln <- VlnPlot(
    seurat_data,
    features.plot = c("nGene", "nUMI", "percent_mito"),
    nCol = 3
)
ggsave(
    filename = "../data/figures/UMI_Violin.pdf",
    plot = umi_vln,
    width = 16, height = 9, device = "pdf"
)

# WT filter

GenePlot(
    seurat_data,
    gene1 = "nUMI",
    gene2 = "nGene",
    cell.ids = grep("WT", colnames(seurat_data@raw.data), value = TRUE)
)
GenePlot(
    seurat_data,
    gene1 = "nUMI",
    gene2 = "percent_mito",
    cell.ids = grep("WT", colnames(seurat_data@raw.data), value = TRUE)
)

seurat_filter_wt <- FilterCells(
    seurat_data,
    subset.names = c("nUMI", "percent_mito"),
    low.thresholds = c(50, -Inf),
    high.threshold = c(1400, 0.1),
    cells.use = WhichCells(seurat_data, ident = "WT")
)

# HET filter

GenePlot(
    seurat_data,
    gene1 = "nUMI",
    gene2 = "nGene",
    cell.ids = grep("HET", colnames(seurat_data@raw.data), value = TRUE)
)
GenePlot(
    seurat_data,
    gene1 = "nUMI",
    gene2 = "percent_mito",
    cell.ids = grep("HET", colnames(seurat_data@raw.data), value = TRUE)
)

seurat_filter_het <- FilterCells(
    seurat_data,
    subset.names = c("nUMI", "percent_mito"),
    low.thresholds = c(0, -Inf),
    high.threshold = c(1250, 0.15),
    cells.use = WhichCells(seurat_data, ident = "HET")
)

# Normalize data, scale, and find variable genes

scale_wt <- NormalizeData(seurat_filter_wt) %>%
    ScaleData(., vars.to.regress = c("nUMI", "percent_mito")) %>%
    FindVariableGenes(., do.plot = FALSE)

scale_het <- NormalizeData(seurat_filter_het) %>%
    ScaleData(., vars.to.regress = c("nUMI", "percent_mito")) %>%
    FindVariableGenes(., do.plot = FALSE)

# Gene selection for CCA (intersect top 1000 per dataset present in both)

genes_use <- unique(c(
    head(rownames(scale_wt@hvg.info), 1000), 
    head(rownames(scale_het@hvg.info), 1000)
)) %>%
    intersect(., rownames(scale_wt@scale.data)) %>%
    intersect(., rownames(scale_het@scale.data))

# CCA

cul3_integrated <- RunCCA(
    scale_wt, scale_het, genes.use = genes_use, num.cc = 30
)
saveRDS(cul3_integrated, "../data/initialCCA.rds")

# Visualze CCA

cca_p1 <- DimPlot(
    object = cul3_integrated, 
    reduction.use = "cca", group.by = "orig.ident", 
    pt.size = 0.5, do.return = TRUE
)
cca_p2 <- VlnPlot(
    object = cul3_integrated, 
    features.plot = "CC1", group.by = "orig.ident", 
    do.return = TRUE
)
ggsave(
    filename = "../data/figures/cca_plots.pdf",
    plot = plot_grid(cca_p1, cca_p2),
    width = 16, height = 9, device = "pdf"
)

bicor_pt <- MetageneBicorPlot(
    cul3_integrated, grouping.var = "orig.ident", 
    dims.eval = 1:30, display.progress = TRUE
)
ggsave(
    filename = "../data/figures/bicor_plot.pdf",
    plot = bicor_pt,
    width = 16, height = 9, device = "pdf"
)

dim_heat <- DimHeatmap(
    object = cul3_integrated, reduction.type = "cca", 
    cells.use = 500, dim.use = 1:9, do.balanced = TRUE, do.return = TRUE
)

# Align subspaces

cul3_aligned <- AlignSubspace(
    cul3_integrated, reduction.type = "cca", 
    grouping.var = "orig.ident", dims.align = 1:20
)
saveRDS(cul3_aligned, "../data/ccaAlign.rds")

# Integrated Analysis

# TSNE

cul3_tsne <- RunTSNE(
    object = cul3_aligned, reduction.use = "cca.aligned", dims.use = 1:20
) %>%
    FindClusters(
        ., reduction.type = "cca.aligned", resolution = 1.5, dims.use = 1:20
    )
tsne_int_plt <- plot_grid(
    TSNEPlot(cul3_tsne, do.return = TRUE, pt.size = 0.5, group.by = "orig.ident"),
    TSNEPlot(cul3_tsne, do.label = T, do.return = T, pt.size = 0.5)
)
saveRDS(cul3_tsne, "../data/cul3TSNE.rds")

# Get markers

cluster_markers <- lapply(
    levels(cul3_tsne@ident),
    function(i) {
        FindConservedMarkers(
            cul3_tsne, ident.1 = i, grouping.var = "orig.ident", print.bar = FALSE
        ) %>%
            mutate(
                Gene = rownames(.)
            )
    }
) %>%
    setNames(., levels(cul3_tsne@ident))
writexl::write_xlsx(cluster_markers, "../data/markers.xlsx")

# Cell deconvolution

cluster_ids <- c(
    "0" = "EN1",
    "1" = "Unknown1",
    "2" = "EN2",
    "3" = "Oligodendrocytes", 
    "4" = "EN3", 
    "5" = "HKG",
    "6" = "EN4",
    "7" = "EN + IN",
    "8" = "Unknown2",
    "9" = "Unknown3",
    "10" = "Astrocytes",
    "11" = "EN5",
    "12" = "Vip+ IN",
    "13" = "IN",
    "14" = "Microglia",
    "15" = "EN6 (?)",
    "16" = "Pericytes",
    "17" = "Unknown4"
)

# Individual genotypes

sc_wt <- StashIdent(cul3_tsne, save.name = "clusterID") %>%
    SetAllIdent(., id = "orig.ident") %>%
    SubsetData(., ident.use = "WT") %>%
    SetAllIdent(., id = "clusterID")
sc_het <- StashIdent(cul3_tsne, save.name = "clusterID") %>%
    SetAllIdent(., id = "orig.ident") %>%
    SubsetData(., ident.use = "HET") %>%
    SetAllIdent(., id = "clusterID")

tsne_plots <- list(
    MERGE = TSNEPlot(cul3_tsne, do.label = T, do.return = T, pt.size = 0.5) +
        ggtitle("MERGE") +
        theme(legend.position = "none"),
    WT = TSNEPlot(sc_wt, do.label = T, do.return = T, pt.size = 0.5) +
        ggtitle("WT") +
        theme(legend.position = "none"),
    HET = TSNEPlot(sc_het, do.label = T, do.return = T, pt.size = 0.5) +
        ggtitle(bquote("Cul3"^"+/-")) +
        theme(legend.position = "none")
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

tsne_plots <- list(
    MERGE = tsne_gg(cul3_tsne, cluster_labels = cluster_ids) +
        ggtitle("MERGE") +
        theme(legend.position = "none"),
    WT = tsne_gg(sc_wt, cluster_labels = cluster_ids) +
        ggtitle("WT") +
        theme(legend.position = "none"),
    HET = tsne_gg(sc_het, cluster_labels = cluster_ids) +
        ggtitle(bquote(bold("Cul3"^"+/-"))) +
        theme(legend.position = "none"),
    LEGEND_IDS = get_legend(
        tsne_gg(cul3_tsne) +
            scale_colour_manual(values = rep("white", length(cluster_ids))) + 
            guides(colour = guide_legend(override.aes = list(size = 3))) +
            theme(
                legend.text.align = 1, 
                legend.text = element_text(size = 18),
                legend.title = element_blank()
            )
        ),
    LEGEND_LABS = get_legend(
        tsne_gg(cul3_tsne, cluster_ids) +
            guides(colour = guide_legend(override.aes = list(size = 3))) +
            theme(
                legend.text = element_text(size = 18), 
                legend.title = element_blank()
            )
        )
)

ggsave(
    filename = "../data/figures/tsne.pdf",
    plot = plot_grid(
        tsne_plots[["MERGE"]], tsne_plots[["WT"]],
        tsne_plots[["HET"]], plot_grid(
            ggplot(),
            tsne_plots[["LEGEND_IDS"]], tsne_plots[["LEGEND_LABS"]],
            ggplot(),
            ncol = 4, nrow = 1, rel_widths = c(0.3, 0.1, 0.1, 0.2)
        ), 
        nrow = 2, ncol = 2
    ),
    width = 16, height = 9, device = "pdf"
)

# Relative Cell Counts

cul3_labeled <- StashIdent(cul3_tsne, save.name = "clusterID") %>%
    AddMetaData(
        ., 
        metadata = data.frame(
            clusterLabel = cluster_ids[match(as.character(.@ident), names(cluster_ids))],
            cell = rownames(.@meta.data),
            row.names = "cell"
        ),
        col.name = "clusterLabel"
    )
total_cell_counts <- setNames(
    sapply(
        levels(cul3_labeled@meta.data$orig.ident),
        function(gn) {
            nrow(filter(cul3_labeled@meta.data, orig.ident == gn))
        }
    ),
    levels(cul3_labeled@meta.data$orig.ident)
)
ggplot(
    data = cul3_labeled@meta.data,
    aes(
        x = reorder(clusterLabel, as.numeric(clusterID)),
        group = orig.ident, fill = orig.ident
    )
) +
    geom_bar(aes(y = ..prop..), colour = "black", position = "dodge") +
    theme(
        axis.text.x = element_text(angle = -30, hjust = 0),
        axis.title = element_blank()
    )

# Neuronal Markers Odds Ratio Barchart

data <- as.matrix(cul3_tsne@data)
meta <- cul3_tsne@meta.data

Inhibitory <- c(
    "Pvalb","Sst","Vip","Htr3a","Ndnf","Calb2","Nos1",
    "Chrna2","Nkx2-1","Gad2","Reln"
)
Excitatory <- c(
    "Camk2a","Cux2","Nr5a1","Scnn1a","Rorb","Rbp4","Ntsr1","Ctgf","Slc17a6",
    "Foxp2","Tmem215","Bend5","Gprin3"
)
markers <- c(Inhibitory, Excitatory)
markers.expr <- data.frame(row.names = colnames(data))[1:ncol(data),]

for (i in markers){
    tmp <- data.frame(
        tryCatch(
            data[i, ],
            error = function(e) rep(0, ncol(data))
        ),
        row.names = colnames(data)
    )
    colnames(tmp) <- i
    markers.expr <- cbind(markers.expr,tmp)
}

table(rownames(meta) == rownames(markers.expr))
markers.expr$genotype <- meta$orig.ident[match(rownames(markers.expr), rownames(meta))]

result <- data.frame()
for (i in 1:(ncol(markers.expr) - 1)) {
    q <- nrow(markers.expr[markers.expr$genotype=="WT" & markers.expr[,i]>0,])
    m <- nrow(markers.expr[markers.expr$genotype=="WT" & markers.expr[,i]==0,])
    k <- nrow(markers.expr[markers.expr$genotype=="HET" & markers.expr[,i]>0,])
    t <- nrow(markers.expr[markers.expr$genotype=="HET" & markers.expr[,i]==0,])
    print(sum(q,m,k,t))
    fisher.out <- fisher.test(matrix(c(k, t, q, m), 2, 2),conf.int=TRUE)
    OR <- fisher.out$estimate
    pval <- fisher.out$p.value
    upCI <- fisher.out$conf.int[1]
    downCI <- fisher.out$conf.int[2]
    result <- rbind(
        result,
        data.frame(
            marker = colnames(markers.expr)[i],
            wt.percent = q / (q + m),
            het.percent = k / (k + t),
            OR = OR,
            upCI = upCI,
            downCI = downCI,
            pval = pval
        )
    )
}

result$type <- rep(
    c("Inhibitory", "Excitatory"),
    c(length(Inhibitory), length(Excitatory))
)
result <- result[result$wt.percent > 0 | result$het.percent > 0, ]
result$fdr <- p.adjust(result$pval, method = "fdr")
result$marker <- factor(result$marker, levels = unique(as.character(result$marker)))

result <- result %>%
    arrange(type, OR) %>%
    mutate(Order = seq(1, nrow(.))) %>%
    filter(type != "NonNeuronal") %>%
    filter(marker != "Rorb") %>%
    mutate(marker = as.character(marker)) %>%
    mutate(SigText = sapply(fdr, function(p) {
        if (p <= 0.001) { return("***") }
        else if (p <= 0.01) { return("**") }
        else if (p <= 0.05) { return("*") }
        else { return("") }
    }))

marker_barplot <- ggplot(
    data = result,
    aes(
        x = reorder(marker, -Order),
        y = OR,
        fill = type
    )
) +
    geom_bar(stat = "identity",width = 0.7) +
    geom_hline(yintercept = 1,linetype = 2, size=1) +
    geom_text(aes(label = SigText), hjust = 0, vjust = 0.7, nudge_y = 0.005, size = 14) +
    scale_y_continuous(
        breaks = c(0, 1, 2),
        expand = expand_scale(mult = c(0, 0.15))
    ) +
    labs(
        x = "",
        y = bquote("Odds Ratio (Cul3"^"+/-"~"vs WT)")
    ) +
    scale_fill_manual(
        values = c(
            "Excitatory" = "firebrick1",
            "Inhibitory" = "springgreen4",
            "***: FDR <= 0.001" = "white",
            "**:  FDR <= 0.01" = "white",
            "*:   FDR <= 0.05" = "white"
        ),
        limits = c(
            "Excitatory",
            "Inhibitory",
            "***: FDR <= 0.001",
            " **: FDR <= 0.01",
            "  *: FDR <= 0.05"
        )
    ) +
    guides(fill = guide_legend(title = "", label.position = "left")) +
    theme(
        text = element_text(size = 22, face = "bold"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, size = 30, face = "bold"),
        axis.text.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(size = 40, face = "bold"),
        legend.justification = c(1, 1),
        legend.position = c(1, 1),
        legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 30)
    ) +
    coord_flip()
marker_barplot

ggsave(
    filename = "../data/figures/OR_neuronal_markers.pdf",
    plot = marker_barplot,
    width = 16, height = 9, device = "pdf"
)

