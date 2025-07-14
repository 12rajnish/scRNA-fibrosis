library(Seurat)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
mito_gene_file <- args[1]

if (file.exists(mito_gene_file)) {
  rabbit_mitochondrial_genes <- readLines(mito_gene_file)
  if (length(rabbit_mitochondrial_genes) == 0) {
    stop("Error: No mitochondrial genes found. Exiting script.")
  }
  print(paste("Loaded", length(rabbit_mitochondrial_genes), "mitochondrial genes from annotation"))
} else {
  stop("Error: Mitochondrial gene list file not found. Exiting script.")
}

data_dir <- "path_to_data/"
output_dir <- "path_to_output/"
dir.create(output_dir, showWarnings = FALSE)

sample_list <- c("AB-1", "AB-2", "AB-3", "Naive-1", "Naive-2", "Naive-3")

for (sample in sample_list) {
  print(paste("Processing:", sample))
  
  sample_path <- file.path(data_dir, sample, "outs", "filtered_feature_bc_matrix")
  if (!dir.exists(sample_path)) {
    print(paste("Skipping:", sample, "- Directory not found"))
    next
  }
  
  seurat_obj <- Read10X(data.dir = sample_path)
  seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = sample)
  
  mito.genes <- rownames(seurat_obj)[rownames(seurat_obj) %in% rabbit_mitochondrial_genes]
  
  if (length(mito.genes) > 0) {
    print(paste("Found", length(mito.genes), "mitochondrial genes for", sample))
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, features = mito.genes)
  } else {
    print(paste("No mitochondrial genes found for", sample))
    seurat_obj[["percent.mt"]] <- rep(0, ncol(seurat_obj))
  }
  
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- SetAssayData(seurat_obj, layer = "data", new.data = GetAssayData(seurat_obj, layer = "counts"))
  
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
  
  saveRDS(seurat_obj, file = file.path(output_dir, paste0(sample, "_filtered.rds")))
  
  print(paste("Completed Processing:", sample))
}
