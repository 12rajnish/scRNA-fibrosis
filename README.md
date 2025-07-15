# scRNA-Alkali-Burn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Author: Rajnish Kumar, Nishant R. Sinha & Rajiv R. Mohan
Copyright(C)2025-2026. Rajnish Kumar, Nishant R. Sinha & Rajiv R. Mohan   
University of Missouri, Columbia          
All Rights Reserved.	             

%  Created by
Rajnish Kumar & Rajiv R. Mohan
Harry S. Truman Memorial Veterans' Hospital, Columbia, Missouri, United States.
Department of Ophthalmology, Veterinary Medicine & Surgery, College of Veterinary Medicine,
University of Missouri, Columbia, Missouri, United States.
Department of Ophthalmology, School of Medicine, University of Missouri, Columbia, Missouri, United States.

%  For more information, contact:

%    Dr. Rajnish Kumar
     University of Missouri-Columbia
     Columbia, MO 65211
     rajnish.kumar@missouri.edu
or
     Dr. Rajiv R. Mohan
     Curators’ Distinguished Professor Ophthalmology and Molecular Medicine
     School of Medicine & College of Veterinary Medicine
     University of Missouri, Columbia
     Columbia, MO 65211
     mohanr@missouri.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scRNA-seq Processing Pipeline (10x Genomics + Cell Ranger + Seurat)

This README provides step-by-step instructions for processing rabbit single-cell RNA-seq data using Cell Ranger and Seurat in R.

Step 1: Install Cell Ranger
Download and extract Cell Ranger:

wget -O cellranger-9.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.0.tar.gz?[TOKEN]"
tar -xzf cellranger-9.0.0.tar.gz

Add it to your PATH or use `which cellranger` to locate.

Step 2: Prepare Reference Genome

(A) Filter for protein-coding genes:
cellranger mkgtf path_to/genomic.gtf path_to/genomic_file.gtf --attribute=gene_biotype:protein_coding

(B) Build genome reference:
cellranger mkref --genome=Rabbit_NZW \
  -fasta=path_to/genomic.fna \
  -genes=path_to/genomic.filtered.gtf

Step 3: Prepare FASTQ Files

If required, rename FASTQ files in 10x format:
e.g.
mv AB-1_S10_R1_merged.fastq.gz AB-1_S10_S1_L001_R1_001.fastq.gz
mv AB-1_S10_R2_merged.fastq.gz AB-1_S10_S1_L001_R2_001.fastq.gz

Step 4: Count Using Cell Ranger

cellranger count --id=All_Samples \
  -transcriptome=path_to/Rabbit_NZW \
  -fastqs=path_to/fastq_files \
  -sample=AB-1,AB-2,AB-3,Naive-1,Naive-2,Naive-3 \
  -create-bam=true

Step 5: Mitochondrial Gene Filtering

Identify mitochondrial genes:

grep -i "mitochondrial" path_to/genomic.gtf | awk -F 'gene_id "' '{print $2}' | awk -F '"' '{print $1}' | sort | uniq > mitochondrial_genes.txt

Step 6: Downstream Processing Scripts

Run the following:

chmod +x process_scRNAseq.sh
./process_scRNAseq.sh

Step 7: Run in R:

`process_scRNAseq.R` for filtering of healthy cells


Step 8: Seurat scRNA-seq Sample Merging and processing Scripts

1. Sample Merging Script (`merge_samples.R`)

This script loads and merges individual Seurat RDS objects representing different biological samples (e.g., Naive and AB groups) into combined Seurat objects.
Main Steps
1. Define file paths for Naive and AB samples
2. Load each sample from disk
3. Convert assays to RNA assay format if needed
4. Merge samples into one Seurat object per group (Naive and AB)
5. Save the merged Seurat objects
Input
- Individual Seurat RDS files located in `path_to_filtered_rds_files/`, such as:
  - Naive-1.rds
  - Naive-2.rds
  - Naive-3.rds
  - AB-1.rds
  - AB-2.rds
  - AB-3.rds
Output: 'All_Samples_merged.rds'
Required Packages
Seurat
ggplot2
dplyr
patchwork
future

Step 9: Run `seurat_code.R` for normalization, variable features, scaling, UMAP, and DEG analysis

This script performs standard preprocessing of a merged single-cell RNA-seq dataset using Seurat.

1. Load merged Seurat object
2. Normalize, identify variable features, and scale
3. Run PCA, clustering, and UMAP (if not already present)
4. Generate a cell count summary per sample and cluster
5. Save UMAP plot with cluster labels
6. Join assay layers and identify cluster-specific marker genes (DEGs)
7. Save processed Seurat object and marker genes table

Input
- A merged Seurat RDS file (`All_Samples_merged.rds`) 

Outputs
Clustered Seurat object: `All_Samples_Clusters.rds`
Cell count summary: `All_Samples_Cluster_CellCounts.csv`
UMAP plot: `All_Samples_UMAP.png`
DEG marker table: `DEGs_Markers.csv`

Required Packages
Seurat
dplyr
tidyr
ggplot2

Step 10: Subclustering Analysis
Run 'subclustering.R'

Key Steps:
1. Clusters subset Selection
2. Preprocessing
3. Dimensionality Reduction & Clustering
4. Marker Gene Detection
5. Visualization (UMAP and feature plots for top DEGs are saved)
6. Subcluster Annotation
   

Outputs
UMAP plot
DEG tables
updated Seurat.rds 

Required Packages
Seurat
ggplot2
dplyr
patchwork

Step 11: Pseudotime and Trajectory Analysis (Monocle3)

Run 'trajectory_and_pseudotime.R'
Steps:
1. Load Seurat object and subset specified clusters.
2. Convert to Monocle3 CellDataSet and inject UMAP coordinates.
3. Learn principal graph and order cells in pseudotime.
4. Generate and save trajectory and pseudotime plots.

Input
Subcluster_Seurat.rds (Seurat object with clusters and UMAP)

Output
Trajectory_SelectedClusters.tiff
Trajectory_Pseudotime_SelectedClusters.tiff

Required Packages
Seurat
monocle3
SeuratWrappers
