suppressPackageStartupMessages({
    library(fgsea)
    library(tibble)
    library(reshape2)
    library(readxl)
    library(org.Hs.eg.db)
    library(msigdbr)
    library(dplyr)
    library(ggplot2)
    library(clusterProfiler)
    library(enrichplot)
    library(dplyr)
    library(DOSE)
    library(tibble)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
})

load_msigdbr <- function() {
    msigdbr_C2 <- msigdbr(species = "human", category = "C2", subcategory="CGP") %>% dplyr::select(gs_name, entrez_gene)
    pathways_C2 = split(x = msigdbr_C2$entrez_gene, f = msigdbr_C2$gs_name)
    msigdbr_H <- msigdbr(species = "human", category = "H") %>% dplyr::select(gs_name, entrez_gene)
    pathways_H = split(x = msigdbr_H$entrez_gene, f = msigdbr_H$gs_name)

    return(list("msigdbr_C2"=msigdbr_C2, "msigdbr_H"=msigdbr_H, "pathways_C2"=pathways_C2, "pathways_H"=pathways_H))
}

run_go_analysis <- function(gene_names, pAdjustMethod = "none") {
    ego <- enrichGO(gene = gene_names, 
                    OrgDb = org.Hs.eg.db, 
                    ont = "CC", 
                    pAdjustMethod = pAdjustMethod, 
                    pvalueCutoff  = 0.05, 
                    readable = TRUE)
    
    p <- dotplot(ego, showCategory=30) + ggtitle("dotplot for GO")

    return(list("ego"=ego, "plot"=p))
} 


run_gseago <- function(ranked_genes, pathway, ont="CC") {
    gsea_res <- gseGO(
        geneList = ranked_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = ont,
        pAdjustMethod = "BH", 
        pvalueCutoff  = 0.05, 
        minGSSize = 10,
        maxGSSize = 300,
        nPermSimple = 100000,
        eps = 1e-20,
        by = "fgsea",
        verbose = FALSE
    )

    return(gsea_res)
}


run_fgsea <- function(RANKS_data, pathway_data) {
    options(repr.plot.height=3, repr.plot.width=10)

    GSEAres <- fgsea(pathways = pathway_data, stats = RANKS_data, minSize = 10, maxSize = 500, nproc = 1)
    GSEAres <- GSEAres %>% arrange(pval)
    GSEAres_subset <- subset(GSEAres, subset=(pval<=0.05))

    p <- ggplot(GSEAres_subset, aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill = NES > 0)) +
        scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
        coord_flip() +
        labs(x = "Pathway", y = "Normalized Enrichment Score",
            title = "Hallmark pathways NES from GSEA") + 
            theme_minimal() +
            theme(legend.position = "none")
    return(p)
}


run_enrichment <- function(marker_results, pathway_data, ont="CC") {
    marker_list <- split(marker_results, marker_results$cluster)
    go_output = list()
    for (cluster in names(marker_list)) {
        df <- marker_list[[cluster]]
        gene_list <- df %>%
            dplyr::filter(!is.na(avg_log2FC)) %>%
            dplyr::filter(!is.na(p_val_adj)) %>%
            dplyr::filter(p_val_adj <= 0.05) %>%
            arrange(desc(avg_log2FC)) %>%
            distinct(gene, .keep_all = TRUE)
        
        gene_list$entrez <- mapIds(
            org.Hs.eg.db,
            keys = gene_list$gene,
            column = "ENTREZID",
            keytype = "SYMBOL",
            multiVals = "first"
        )
        
        ranked_genes <- gene_list$avg_log2FC
        names(ranked_genes) <- gene_list$entrez
        ranked_genes <- sort(ranked_genes, decreasing = TRUE)
        ranked_genes <- ranked_genes[!is.na(names(ranked_genes))]
        
        out <- run_gseago(ranked_genes, pathway_data, ont)
        print(paste0("Cluster ", as.character(cluster)))
        go_output[[cluster]] <- out
    }
    
    return(go_output)
}




plot_gene_heatmap <- function(seurat_obj, marker_list, cell_order=NULL, break_lim=2) {

    avg_exp_mat <- AverageExpression(
        seurat_obj, 
        assays = "RNA", 
        slot = "data"
    )$RNA
    
    celltype_cols <- names(marker_list)
    if (!is.null(cell_order)) {
        cell_order = colnames(avg_exp_mat)
    }
    
    avg_exp_mat <- avg_exp_mat[, cell_order]
    
    avg_exp_df <- avg_exp_mat %>%
        as.data.frame() %>%
        rownames_to_column("gene")

    all_markers_in_order <- unlist(marker_list, use.names = FALSE)
    all_markers_in_order <- unique(all_markers_in_order)

    plot_df <- avg_exp_df %>% 
        dplyr::filter(gene %in% all_markers_in_order) %>%
        slice(match(all_markers_in_order, gene))

    plot_mat <- as.matrix(plot_df[, cell_order])
    
    rownames(plot_mat) <- plot_df$gene
    
    marker_counts <- table(factor(unlist(marker_list), levels = rownames(plot_mat)))
    marker_counts <- marker_counts[marker_counts > 0]
    gaps_row_positions <- cumsum(as.numeric(marker_counts))

    plot <- pheatmap(
        plot_mat,
        scale="row",
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        fontsize_col = 12,
        fontsize_row = 12,
        border = TRUE,
        color = colorRampPalette(c("dodgerblue4","white","darkred"))(100),
        border_color = "white",
        show_rownames = TRUE,
        show_colnames = TRUE,
        main = "Marker Heatmap"
    )

    return(plot)
    
}


assign_dominant_cluster <- function(row, plot_mat) {
    cell_type <- colnames(plot_mat)[which.max(row)]
    return(cell_type)
}




plot_pathway_heatmap <- function(gsea_output, cell_order, left_margin=unit(2, "cm")) {
    plot_data_list <- list()
    for (cluster in names(gsea_output)) {
        res <- as.data.frame(gsea_output[[cluster]])
        if (nrow(res) > 0) {
            res <- res %>% 
                dplyr::filter(NES != 0) %>% 
                dplyr::filter(p.adjust <= 0.05)
            res$Cluster <- cluster
            plot_data_list[[cluster]] <- res[, c("Description", "NES", "Cluster")]
        }
    }
    
    plot_data_combined <- do.call(rbind, plot_data_list)
    plot_mat <- reshape2::dcast(plot_data_combined, Description ~ Cluster, value.var = "NES")
    
    rownames(plot_mat) <- plot_mat$Description
    plot_mat$Description <- NULL
    plot_mat[is.na(plot_mat)] <- 0
    
    valid_cell_order <- intersect(cell_order, colnames(plot_mat))
    if (length(valid_cell_order) == 0) {
        stop("None of the provided cell_order values match the columns in the data matrix.")
    }
    plot_mat <- plot_mat[, valid_cell_order, drop = FALSE]
    
    dominant_cluster <- apply(plot_mat, 1, assign_dominant_cluster, plot_mat = plot_mat)
    
    pathway_order_df <- data.frame(
        Pathway = rownames(plot_mat),
        DominantCellType = dominant_cluster,
        MaxNES = apply(plot_mat, 1, max),
        stringsAsFactors = FALSE
    )
    
    pathway_order_df$DominantCellType <- factor(pathway_order_df$DominantCellType, levels = valid_cell_order)
    pathway_order_df <- pathway_order_df %>% arrange(DominantCellType, desc(MaxNES))
    plot_mat_ordered <- plot_mat[pathway_order_df$Pathway, , drop = FALSE]
    plot_mat_ordered <- plot_mat_ordered[rowSums(plot_mat_ordered > 0) > 0, , drop = FALSE]
    
    plot <- pheatmap(
        plot_mat_ordered,
        scale = "row",
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        fontsize_col = 24,
        fontsize_row = 12,
        cell_height = 10,
        border = TRUE,
        lhei = c(1, 10), 
        lwid = c(1, 10),
        color = colorRampPalette(c("dodgerblue4", "white", "darkred"))(100),
        border_color = "white",
        show_rownames = TRUE,
        show_colnames = TRUE,
        fontfamily = "sans",
        angle_col = 90,
        silent = TRUE,
        main = ""
    )
    


    plot$gtable$widths[2] <- plot$gtable$widths[2] + left_margin


    return(plot)
}
