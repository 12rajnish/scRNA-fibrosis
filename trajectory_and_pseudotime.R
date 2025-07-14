library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
library(patchwork)

seurat_obj <- readRDS("Subcluster_Seurat.rds")

selected_clusters <- c("0", "1", "2", "3")
seurat_obj$seurat_clusters <- as.character(seurat_obj$seurat_clusters)
seurat_obj <- subset(seurat_obj, idents = selected_clusters)

cds <- as.cell_data_set(seurat_obj)
rowData(cds)$gene_short_name <- rownames(cds)
reducedDims(cds)$UMAP <- seurat_obj@reductions$umap@cell.embeddings
cds@colData$seurat_clusters <- seurat_obj$seurat_clusters
cds@clusters$UMAP$clusters <- factor(seurat_obj$seurat_clusters)
cds@clusters$UMAP$partitions <- setNames(factor(rep(1, ncol(cds))), colnames(cds))

cds <- preprocess_cds(cds, num_dim = 50)
cds <- learn_graph(cds)
cds <- order_cells(cds)

dir.create("Trajectory_Plots_SelectedClusters", showWarnings = FALSE, recursive = TRUE)

tiff("Trajectory_SelectedClusters.tiff", width = 10, height = 8, units = "in", res = 300)
plot_cells(cds, color_cells_by = "seurat_clusters", label_groups_by_cluster = TRUE,
           label_leaves = FALSE, label_branch_points = TRUE)
dev.off()

tiff("Trajectory_Pseudotime_SelectedClusters.tiff", width = 10, height = 8, units = "in", res = 300)
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE,
           label_leaves = FALSE, label_branch_points = TRUE)
dev.off()

set.seed(42)
cell_ids <- sample(colnames(cds), size = min(3000, ncol(cds)))
cds_sub <- cds[, cell_ids]

tiff("Trajectory_Pseudotime_Downsampled_SelectedClusters.tiff", width = 10, height = 8, units = "in", res = 300)
plot_cells(cds_sub, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE,
           label_branch_points = TRUE, graph_label_size = 6, cell_size = 0.5)
dev.off()

genes_to_plot <- c("TGFB1", "PDGF", "VEGFA", "FGF2", "COL1A1", "JAM3", "LAMA2")
genes_found <- genes_to_plot[genes_to_plot %in% rownames(cds)]

if (length(genes_found) > 0) {
  tiff("Gene_Trend_SelectedClusters.tiff", width = 10, height = 6, units = "in", res = 300)
  plot_genes_in_pseudotime(cds[genes_found, ], color_cells_by = "pseudotime")
  dev.off()
}

deg_pseudo <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
significant_genes <- subset(deg_pseudo, q_value < 0.05)

modules <- find_gene_modules(cds[rownames(significant_genes), ], resolution = 1e-2)

tiff("Gene_Modules_SelectedClusters.tiff", width = 12, height = 8, units = "in", res = 300)
plot_cells(cds, genes = modules, label_cell_groups = FALSE, show_trajectory_graph = FALSE)
dev.off()

