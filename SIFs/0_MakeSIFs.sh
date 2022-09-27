#!/bin/bash -l

# Record of how to pull and build singularity images for modules that have gone out of date or are too hard to reinstall

module load singularity/3.9.3

# Trim Galore Image
singularity build trim-galore_0.6.5.sif docker://quay.io/biocontainers/trim-galore:0.6.5--0

# RSeQC
singularity build RSeQC_4.0.0.sif docker://quay.io/biocontainers/rseqc:4.0.0--py310h1425a21_2

# MultiQC Image
singularity build multiQC_1.13.sif docker://quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0
