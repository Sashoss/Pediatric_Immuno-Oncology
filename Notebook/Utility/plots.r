library(ggplot2)


create_multiple_dotplots <- function(seurat_data, marker_list) {
    dotplot_list <- list()
    for (cell_name in names(marker_list)) {
        dotplot <- DotPlot(seurat_data, features = marker_list[[cell_name]]) +
            RotatedAxis() +
            scale_colour_gradientn(colours = c("blue", "white", "red"), space = "Lab") +
            theme_minimal() +
            theme(
                legend.position = "right",
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
                axis.text.y = element_text(size = 14)
            ) +
            ggtitle(paste(cell_name))
    
        dotplot_list[[cell_name]] <- dotplot
    }
    return(dotplot_list)
}


create_single_dotplot <- function(seurat_obj, scores_names, group_by) {
    if (is.null(seurat_obj)) {
        stop("Error: The provided Seurat object is NULL.")
    }
    dotplot <- DotPlot(
        seurat_obj,
        features = scores_names,
        group.by = group_by
    ) +
    scale_color_gradientn(colors = c("blue", "white", "red")) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    labs(
        x = "Cell Type Signature",
        y = "Cluster",
        color = "Module Score"
    )

    return(dotplot)
}
