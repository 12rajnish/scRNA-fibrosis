library(Seurat)

all_samples <- c("AB-1", "AB-2", "AB-3", 
                 "Naive-1", "Naive-2", "Naive-3")

merge_samples <- function(sample_list, merged_name) {
  merged_seurat <- NULL

  for (sample_name in sample_list) {
    print(paste("Loading:", sample_name))

    filtered_data_path <- paste0("path_to_filtered_rds/", sample_name, "_.rds")

    if (!file.exists(filtered_data_path)) {
      print(paste("File not found for", sample_name, "- Skipping"))
      next
    }

    seurat_obj <- readRDS(filtered_data_path)
    seurat_obj$sample <- sample_name
    seurat_obj <- RenameCells(seurat_obj, new.names = paste0(sample_name, "_", colnames(seurat_obj)))

    if (is.null(merged_seurat)) {
      merged_seurat <- seurat_obj
    } else {
      merged_seurat <- merge(merged_seurat, y = seurat_obj)
    }
  }

  if (!is.null(merged_seurat)) {
    merged_data_path <- paste0("path_to_filtered_rds/", merged_name, "_merged.rds")
    saveRDS(merged_seurat, file = merged_data_path)
    print(paste("Merged dataset saved as:", merged_data_path))
  } else {
    print(paste("No valid samples found for merging:", merged_name))
  }

  return(merged_seurat)
}

all_merged <- merge_samples(all_samples, "All_Samples")
print("Merging completed")
