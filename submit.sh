#!/usr/bin/env bash

#SBATCH -t 20:00:00 
#SBATCH -p shared
#SBATCH --mem=5G
#SBATCH -J SnakeJob
#SBATCH -o mainSnakeJob.out

snakemake -j 999 \
    --cluster-config cluster.yaml \
    --cluster "sbatch \
        -t {cluster.time} \
        -p {cluster.partition} \
        --mem={cluster.mem} \
        -J {cluster.job-name}\
        --output={cluster.output}\
        --error={cluster.error}"
