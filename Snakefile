"""
__author__ = "Kevin Chau"
__date-created__ = "2019-03-02"
__date-modified__ = "2019-03-02"
"""

import os

configfile: "config.yaml"

# Setup pipeline directories
if not os.path.exists("data"):
    os.makedirs("data")
if not os.path.exists("log"):
    os.makedirs("log")

# Declare software
cellranger = config["cellranger"]["path"]

# References
txtome = config["cellranger"]["transcriptome"]

rule all:
   input: 
        directory(expand("data/{sample_id}", sample_id = config["cellranger"]["ids"].keys()))

rule counts:
    input: 
        fqs = lambda wildcard: 
            f"{config['cellranger']['ids'][wildcard.sample_id]['fastq']}"
    output: 
        directory("data/{sample_id}"),
    params:
        samples = lambda wildcard:
            f"{config['cellranger']['ids'][wildcard.sample_id]['samples']}"
    shell: 
        f"{cellranger} count "
        "--id={wildcards.sample_id} "
        "--fastqs={input.fqs} "
        "--sample={params.samples} "
        f"--transcriptome={txtome}; "
        "cp -r {wildcards.sample_id}/outs/filtered_feature_bc_matrix/* data/{wildcards.sample_id}/"

