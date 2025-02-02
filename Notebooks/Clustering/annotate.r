library(readxl)
library(ggplot2)
library(Seurat)



seurat_obj <- readRDS(paste0("./out/integrated_seurat_clusters_0.2.rds"))
print("Combined seurat object")
print(seurat_obj)

seurat_obj <- SetIdent(seurat_obj, value = "RNA_snn_res.0.5")
umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)

ggsave(
  filename = "./out/umap_05.png",
  plot = umap,
  width = 8,
  height = 6
)

markers <- read.table("./in/marker_list.tsv", header=TRUE)
marker_list = list()
for (i in 1:nrow(markers)) {
    cell_type <- markers[i,]$CELL_TYPE
    marker_obj <- strsplit(markers[i,]$MARKERS, split=",")
    marker_list[[cell_type]] <- marker_obj[[1]]
}


create_and_save_dotplot <- function(cell_type, seurat_obj, marker_list, output_dir, umap) {
  options(repr.plot.width = 14, repr.plot.height = 7)
  
  # Check if the cell type exists in the marker list
  if (!cell_type %in% names(marker_list)) {
    stop(paste("Cell type", cell_type, "not found in marker list"))
  }
  
  # Create the DotPlot
  dotplot <- DotPlot(seurat_obj, features = marker_list[[cell_type]]) +
    RotatedAxis() +
    scale_colour_gradientn(colours = c("blue", "white", "red"), space = "Lab") +
    theme_minimal() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)  # Center the title
    ) +
    ggtitle(paste("Expression of Markers Across Cell Types:", cell_type))
  
  combined_plot <- dotplot | umap
  # # Save the plot
  # output_file <- file.path(output_dir, paste0(cell_type, "_dotplot.png"))
  # ggsave(output_file, plot = combined_plot, width = 10, height = 5)
  
  # message(paste("Saved plot for", cell_type, "to", output_file))
return(combined_plot)
}

# Replace "NEURONS" with the desired cell type
cell_types <- c("NEURONS", "ASTROCYTES", "ENDOTHELIAL", "PERICYTES", 
                "MESENCHYMAL", "MICROGLIA", "TAMS", "T", "NK", "B", 
                "TUMOR", "MES-LIKE", "AC-LIKE", "OPC-LIKE", "NPC-LIKE")

output_dir <- "out/dotplots"  # Set the directory to save plots
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

for (cell_type in cell_types) {
  create_and_save_dotplot(cell_type, seurat_obj, marker_list, output_dir, umap)
}

