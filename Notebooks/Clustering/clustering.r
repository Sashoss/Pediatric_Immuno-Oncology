library(Seurat)


seurat_object <- readRDS("../integrate/out/integrated_harmony.rds")

seurat_object <- FindNeighbors(seurat_object, reduction="harmony", dims = 1:30)

resolutions_to_test <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99)
for (res in resolutions_to_test) {
  seurat_object <- FindClusters(seurat_object, resolution = res, verbose = FALSE)
}

library(clustree)

# Plot the cluster tree
p <- clustree(seurat_object, direction = "downwards")
# Save the plot
ggsave(
  filename = "./out/cluster_tree.png",
  plot = p,
  width = 8,
  height = 6
)


library(ggplot2)

Idents(seurat_object) <- seurat_object$RNA_snn_res.0.2

umap_obj <- DimPlot(seurat_object, reduction = "umap") +
  ggtitle("UMAP - Sequencing Type") +
  theme_minimal()

ggsave(umap_obj, filename="./out/cluster_dimplot_0.2.png")

saveRDS(seurat_object, "./out/integrated_seurat_clusters_0.2.rds")