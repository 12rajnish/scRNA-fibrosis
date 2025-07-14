library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

seurat_obj <- readRDS("All_Samples_Clusters.rds")

clusters_to_subcluster <- c(x, y, z) # cluster numbers=x, y, z
sub_cluster_obj <- subset(seurat_obj, idents = clusters_to_subcluster)

sub_cluster_obj <- NormalizeData(sub_cluster_obj)
sub_cluster_obj <- FindVariableFeatures(sub_cluster_obj, selection.method = "vst", nfeatures = 2000)
sub_cluster_obj <- ScaleData(sub_cluster_obj)
sub_cluster_obj <- RunPCA(sub_cluster_obj, npcs = 20)

sub_cluster_obj <- FindNeighbors(sub_cluster_obj, dims = 1:20)
sub_cluster_obj <- FindClusters(sub_cluster_obj, resolution = 0.1)
sub_cluster_obj <- RunUMAP(sub_cluster_obj, dims = 1:20)

umap_plot <- DimPlot(sub_cluster_obj, reduction = "umap", label = TRUE, group.by = "seurat_clusters") +
  ggtitle("Sub-Clustered UMAP")
ggsave("Subcluster_UMAP.tiff", plot = umap_plot, width = 8, height = 6, dpi = 300, device = "tiff")

sub_cluster_markers <- FindAllMarkers(sub_cluster_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sub_cluster_markers, file = "DEGs_Subcluster.csv", row.names = FALSE)

top_genes <- sub_cluster_markers %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC) %>%
  pull(gene)
tiff("FeaturePlot_Subclusters.tiff", width = 10, height = 12, units = "in", res = 300, compression = "lzw")
FeaturePlot(sub_cluster_obj, features = top_genes, cols = c("lightblue", "red4"), reduction = "umap")
dev.off()

seurat_obj$sub_cluster <- as.character(seurat_obj$seurat_clusters)
seurat_obj$sub_cluster[Cells(sub_cluster_obj)] <- paste0("Sub_", sub_cluster_obj$seurat_clusters)

saveRDS(sub_cluster_obj, file = "Subcluster_Seurat.rds")
saveRDS(seurat_obj, file = "Seurat_With_Subclusters.rds")

