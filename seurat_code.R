library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

data_dir <- "path_to_data/"
merged_data_path <- paste0(data_dir, "All_Samples_merged.rds")
clustered_data_path <- paste0(data_dir, "All_Samples_Clusters.rds")
output_csv_path <- paste0(data_dir, "All_Samples_Cluster_CellCounts.csv")
umap_plot_path <- paste0(data_dir, "All_Samples_UMA.png")
degs_csv_path <- paste0(data_dir, "DEGs_Markers.csv")

if (!file.exists(merged_data_path)) {
  stop("Merged dataset file does not exist.")
}
seurat_obj <- readRDS(merged_data_path)

if (!"pca" %in% names(seurat_obj@reductions)) {
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = "percent.mt")
  seurat_obj <- RunPCA(seurat_obj, npcs = 50)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.1)
  saveRDS(seurat_obj, file = clustered_data_path)
}

Idents(seurat_obj) <- seurat_obj$seurat_clusters

cell_counts <- as.data.frame(table(seurat_obj@meta.data$sample, seurat_obj@meta.data$seurat_clusters))
colnames(cell_counts) <- c("Sample_ID", "Cluster", "Freq")
cell_counts_wide <- cell_counts %>%
  pivot_wider(names_from = Cluster, values_from = Freq, values_fill = 0)
cell_counts_wide$Total_Cells <- rowSums(cell_counts_wide[,-1])
cell_counts_wide <- cell_counts_wide %>%
  relocate(Total_Cells, .after = Sample_ID)
write.csv(cell_counts_wide, output_csv_path, row.names = FALSE)

if (!"umap" %in% names(seurat_obj@reductions)) {
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
}

umap_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 5) +
  ggtitle("UMAP Clustering (All Samples)") + theme_minimal()
ggsave(filename = umap_plot_path, plot = umap_plot, width = 8, height = 6)

DefaultAssay(seurat_obj) <- "RNA"
if (!"data" %in% names(seurat_obj[["RNA"]]) || !"counts" %in% names(seurat_obj[["RNA"]])) {
  seurat_obj <- JoinLayers(seurat_obj, layers = c("counts", "data", "scale.data"))
}

if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  stop("'seurat_clusters' column missing in metadata.")
}

Idents(seurat_obj) <- seurat_obj$seurat_clusters
degs_markers <- tryCatch({
  FindAllMarkers(seurat_obj, 
                 only.pos = TRUE, 
                 min.pct = 0.1, 
                 logfc.threshold = 0.25, 
                 test.use = "wilcox")
}, error = function(e) {
  return(NULL)
})

if (!is.null(degs_markers) && nrow(degs_markers) > 0) {
  write.csv(degs_markers, degs_csv_path, row.names = FALSE)
}

saveRDS(seurat_obj, file = clustered_data_path)
