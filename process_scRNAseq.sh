#!/bin/bash

# Define input and output paths
GTF_FILE="path_to/genomic.gtf"
MITO_GENES_FILE="path_to/mitochondrial_genes.txt"
R_SCRIPT="path_to/process_scRNAseq.R"

grep -i "mitochondrial" "$GTF_FILE" | awk -F 'gene_id "' '{print $2}' | awk -F '"' '{print $1}' | sort | uniq > "$MITO_GENES_FILE"

if [ -s "$MITO_GENES_FILE" ]; then
    echo "Mitochondrial gene extraction complete. Found $(wc -l < "$MITO_GENES_FILE") genes."
else
    echo "Warning: No mitochondrial genes found. Check your GTF annotation file."
    exit 1
fi
echo "Running Seurat R script for mitochondrial filtering..."
Rscript "$R_SCRIPT" "$MITO_GENES_FILE"
