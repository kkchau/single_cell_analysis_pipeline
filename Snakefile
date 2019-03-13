"""
__author__ = "Kevin Chau"
__date-created__ = "2019-03-04"
__date-modified__ = "2019-03-09"
"""

import os

localrules: all, seurat_process

configfile: "config.yaml"

# Setup pipeline directories
if not os.path.exists("data"):
    os.makedirs("data")
if not os.path.exists("data/figures"):
    os.makedirs("data/figures")
if not os.path.exists("log"):
    os.makedirs("log")

subworkflow data_setup:
    workdir: "data_setup"
    snakefile: "data_setup/Snakefile"
    configfile: "data_setup/config.yaml"

rule all:
   input: 
        "data/ccaAlign.rds"

rule seurat_process:
    input:
        data_setup(
            "data/{sample_id}/seurat.rds"
        )
    output:
        "data/{sample_id}/seurat_processed.rds"
    params:
        nUMI_low = lambda wildcard:
            f"{config['ids'][wildcard.sample_id]['nUMI']['LOW']}",
        nUMI_high = lambda wildcard:
            f"{config['ids'][wildcard.sample_id]['nUMI']['HIGH']}",
        percent_mito_low = lambda wildcard:
            f"{config['ids'][wildcard.sample_id]['percent_mito']['LOW']}",
        percent_mito_high = lambda wildcard:
            f"{config['ids'][wildcard.sample_id]['percent_mito']['HIGH']}"
    script:
        "scripts/seurat_process.R"

rule seurat_combine:
    input:
        expand("data/{sample_id}/seurat_processed.rds", sample_id = config["ids"].keys())
    output:
        "data/initialCCA.rds"
    script:
        "scripts/seurat_combine.R"

rule seurat_align:
    input:
        "data/initialCCA.rds"
    output:
        alignment_plots = "data/figures/seurat_alignments.pdf",
        cca_align = "data/ccaAlign.rds"
    params:
        heatmap_start = config['heatmap_dimensions']['START'],
        heatmap_end = config['heatmap_dimensions']['END'],
        align_dims_start = config['dim_to_align']['START'],
        align_dims_end = config['dim_to_align']['END']
    script:
        "scripts/seurat_align.R"

rule tsne:
    input:
        "data/ccaAlign.rds"
    output:
        "data/seurat_tsne.rds",
        "data/figures/merged_tsne.pdf",
        expand("data/figures/{sample_id}_tsne.pdf", sample_id = config["ids"].keys())
    params:
        align_dims_start = config['dim_to_align']['START'],
        align_dims_end = config['dim_to_align']['END'],
        tsne_resolution = config['tsne_resolution']