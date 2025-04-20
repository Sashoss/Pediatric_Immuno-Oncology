
library(RColorBrewer)
library(ggplot2)
library(scales)
library(viridis)
library(scales)
library(slingshot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(patchwork)
library(igraph)
library(ggbeeswarm)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(RColorBrewer)

plot_slingshot_trajectory <- function(seurat_obj, slingshot_obj) {

    curves <- slingCurves(slingshot_obj, as.df=T)
    sling_mst_obj <- slingMST(slingshot_obj, as.df=T)
    centers <- as.data.frame(do.call(rbind, V(slingMST(slingshot_obj))$coordinates))
    labels <- V(slingMST(slingshot_obj))$name
    centers$label <- labels

    cell_embeddings <- as.data.frame(Embeddings(seurat_obj, "umap"))
    cell_embeddings$cell_label <- seurat_obj$cell_label 

    cell_plot <- ggplot(cell_embeddings, aes(x = umap_1, y = umap_2, color = cell_label)) +
        geom_point(alpha = 0.6, size = 4)  +  
        theme_classic() +
        ggtitle("Cells with Slingshot MST Overlay")

    cell_plot <- cell_plot +
    geom_path(data = curves %>% arrange(Order),
            aes(x = umap_1, y = umap_2, group = Lineage),
            color = "black", size = 1)

    cell_plot <- cell_plot +
        geom_text(data = centers, 
                  aes(x = umap_1, 
                      y = umap_2, 
                      label = label
                     ), 
                  color = "black", 
                  size = 5
                 ) 

    return(cell_plot)

}



plot_slingshot_lineage <- function(seurat_obj, slingshot_obj) {
    slingshot_obj_test <- slingshot_obj
    colData(slingshot_obj_test)$selected_ident <- as.character(seurat_obj$cell_label)
    slinsghot_data <- colData(slingshot_obj_test)
    column_obj <- grep(pattern="slingPseudotime", x=names(slinsghot_data))
    slingshot_subset <- as.data.frame(slinsghot_data[, column_obj, drop=FALSE])
    slingshot_subset$selected_ident <- slinsghot_data$selected_ident

    slingshot_long <- slingshot_subset %>%
        select(selected_ident, starts_with("slingPseudotime_")) %>% 
        pivot_longer(
            cols = starts_with("slingPseudotime_"),
            names_to = "Lineage",
            values_to = "Pseudotime"
        ) %>%
        drop_na(Pseudotime)

    n_values <- nlevels(as.factor(slingshot_long$selected_ident))
    select_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n_values)

    plot_slingshot_all <- ggplot(slingshot_long, aes(x = Pseudotime, y = Lineage)) +
        geom_quasirandom(
            groupOnX = FALSE,
            size=4,
            aes(color = factor(selected_ident)),
            alpha = 1
        ) +
        scale_color_manual(values = select_colors, name = "Clusters") +
        theme_classic() +
        theme(text = element_text(size = 20)) +
        xlab("Pseudotime") +
        ylab("Lineages") +
        ggtitle("Slingshot Pseudotime by Lineage (Start=MES)")

    return(plot_slingshot_all)

}