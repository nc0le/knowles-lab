#!/bin/bash

#SBATCH -p cpu
#SBATCH --mem=60g
#SBATCH -c 8
#SBATCH --time=40:00:00
#SBATCH --job-name=transcriptome
#SBATCH --output=%A.out
#SBATCH --error=%A.err

nextflow run epi2me-labs/wf-transcriptomes \
        --direct_rna \
        --fastq 'cDNA.fastq' \
        --minimap2_opts '-k 12 -u b -x splice --secondary=no' \
        --ref_genome 'GRCh38.primary_assembly.genome.fa' \
        --stringtie_opts " " \
        --threads 8 \
        -profile singularity
