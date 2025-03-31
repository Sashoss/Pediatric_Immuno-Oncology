# Author - Shiwani Limbu (shiwanilimbu1122@gmail.com, slimbu@kumc.edu)

library(tradeSeq)
library(ggplot2)
library(SingleCellExperiment)
library(writexl)
library(pheatmap)
library(qs)
library(gridExtra)
library(openxlsx)
library(dplyr)
library(msigdbr)
library(fgsea)
library(knitr)
library(clusterProfiler)
library(forcats)
library(cluster)


splitAssocResByCondition <- function(assocRes) {
  global_cols <- c("waldStat", "df", "pvalue", "meanLogFC")
  cols_2p   <- grep("_condition2_percent$",  colnames(assocRes), value = TRUE)
  cols_6p   <- grep("_condition6_percent$",  colnames(assocRes), value = TRUE)
  cols_ctrl <- grep("_conditionCtrl$",       colnames(assocRes), value = TRUE)
  assocRes_2p   <- assocRes[, c(global_cols, cols_2p), drop = FALSE]
  assocRes_6p   <- assocRes[, c(global_cols, cols_6p), drop = FALSE]
  assocRes_ctrl <- assocRes[, c(global_cols, cols_ctrl), drop = FALSE]
  
  retObj <- list(
    `2_percent` = assocRes_2p,
    `6_percent` = assocRes_6p,
    `Ctrl`      = assocRes_ctrl
  )
  return(retObj)
}

substrRight <- function(x, n){
    # pull last n characters of a string
    substr(x, nchar(x)-n+1, nchar(x))
}

get_lineage_genes <- function(assocRes, condition_name, LINEAGE_INDICES) {
    lineage_genes <- list()
    for (lineage_name in names(LINEAGE_INDICES)) {

        lineage_index <- substrRight(lineage_name, 1)
        if (condition_name == "all") {
            pvalue_column_name = paste0("pvalue_", lineage_index)
        } else {
            pvalue_column_name = paste0("pvalue_lineage", lineage_index, "_condition", condition_name)
        }

        lineage_genes[[lineage_name]] <-  rownames(assocRes)[
                    which(p.adjust(assocRes[[pvalue_column_name ]], "fdr") <= 0.05)
        ]
    }

    return(lineage_genes)
}


get_heatSmooth <- function(yhatSmooth, LINEAGE_INDICES) {
    retObj <- list()
    print(dim(yhatSmooth))
    for (name in names(LINEAGE_INDICES)) {
        retObj[[name]] <- pheatmap(
            t(scale(t(yhatSmooth[, LINEAGE_INDICES[[name]]]))),
            cluster_cols = FALSE,
            show_rownames = FALSE, show_colnames = FALSE, main = name,
            legend = FALSE, silent = TRUE
        )
    }

    return(retObj)
}


add_optimal_cluster_to_heatSmooth <- function(yhatSmooth, optimal_k_list, LINEAGE_INDICES) {
    retObj <- list()
    for (name in names(LINEAGE_INDICES)) {
        retObj[[name]] <- pheatmap(
            t(scale(t(yhatSmooth[, LINEAGE_INDICES[[name]]]))),
            cluster_cols = FALSE,
            show_rownames = FALSE, show_colnames = FALSE, main = name,
            cutree_rows = optimal_k_list[[name]],
            legend = FALSE, silent = TRUE
        )
    }
    
    return(retObj)
}





optimal_kmeans <- function(data_matrix, clustering_range = c(2,10), nstart = 25, seed = 123) {
    # Apply silhouette score based calculation of optimal number of kmeans clusters for each lineages
    k_range = clustering_range[1]:clustering_range[2]
    silhouette_scores <- numeric(length(k_range))
    kmeans_models <- vector("list", length(k_range))
    dist_matrix <- stats::dist(data_matrix)
    for (i in seq_along(k_range)) {
        k <- k_range[i]
        set.seed(seed)
        km_res <- kmeans(data_matrix, centers = k, nstart = nstart)
        sil <- cluster::silhouette(km_res$cluster, dist_matrix)
        avg_sil_width <- mean(sil[, 3])
        silhouette_scores[i] <- avg_sil_width
        kmeans_models[[i]] <- km_res
    }

    best_k <- k_range[which.max(silhouette_scores)]
    best_model <- kmeans_models[[which.max(silhouette_scores)]]
    return(best_k)
}



get_lineage_cluster <- function(cluster_genes, 
                                heatSmooth_obj, 
                                yhatSmooth_dataMatrix, 
                                lineage_label, 
                                LINEAGE_INDICES, 
                                clustering_range
                                ) {
    
    dataMatrix_col_index=LINEAGE_INDICES[[lineage_label]]
    data_matrix <- yhatSmooth_dataMatrix[,dataMatrix_col_index[1]:dataMatrix_col_index[2]]
    optimal_k <- optimal_kmeans(data_matrix, clustering_range=clustering_range)
    cl_lineage <- cutree(heatSmooth_obj$tree_row, k = optimal_k)
    for (cluster_id in unique(cl_lineage)) {
        cluster_genes[[lineage_label]][[paste0("cluster_", cluster_id)]] <- names(cl_lineage[cl_lineage == cluster_id])
    }
    return(list("cluster_genes"=cluster_genes,
                "optimal_k"=optimal_k
        )
    )
}


run_go_lineage_clusters <- function(cluster_genes, org_db = "org.Hs.eg.db", p_cutoff = 0.05, q_cutoff = 0.05) {
    go_results <- list()

    for (cluster_name in names(cluster_genes)) {
        genes <- cluster_genes[[cluster_name]]
        if (length(genes) > 0) {
            go_result <- enrichGO(
                gene          = genes,
                OrgDb         = org_db,
                ont           = "BP",
                pAdjustMethod = "BH",
                keyType       = "SYMBOL",
                pvalueCutoff  = p_cutoff,
                qvalueCutoff  = q_cutoff,
                readable      = TRUE
            )
            go_results[[cluster_name]] <- go_result
        }
    }

    return(go_results)
}

run_gsea <- function(assocRes, gene_set, wald_stat_column, m_list, pvalue=0.05) {
    assocRes_data <- assocRes[gene_set,]
    stats <- assocRes[[wald_stat_column]]  
    finite_indices <- is.finite(stats) 
    stats <- stats[finite_indices]      
    names(stats) <- rownames(assocRes)[finite_indices]  
    eaRes <- fgsea(pathways = m_list, stats = stats, minSize = 15, maxSize = 500, scoreType = "pos")
    ooEA <- order(eaRes$pval, decreasing = FALSE)
    fgseaRes <- eaRes[ooEA, ]
    fgseaRes <- fgseaRes[fgseaRes$pval <= pvalue, ]
    fgseaRes <- fgseaRes[!is.na(fgseaRes$pval), ]
    retObj <- list(
        `eaRes` = eaRes,
        `fgseaRes` = fgseaRes
    )
    return(retObj)
}

plot_fgsea <- function(fgsea_result, pvalue=0.05) {

    plot_data <- fgsea_result %>%
        filter(!is.na(pval), pval <= pvalue) %>%
        arrange(desc(NES)) %>%
        slice_head(n = 20) %>%
        mutate(
            pathway          = fct_reorder(pathway, NES),
            minus_log10_pval = -log10(pval)
        )

    p <- ggplot(plot_data, aes(x = minus_log10_pval, y = pathway)) +
        geom_point(aes(size = NES), color = "steelblue") + 
        scale_size_continuous(name = "Normalized Enrichment Score (NES)", range = c(3, 10)) + 
        labs(
            title = "Pathway Enrichment Analysis",
            x = expression(-log[10](p-value)),
            y = "",
            caption = ""
        ) +
        theme_minimal(base_size = 14) + 
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"), 
            axis.text.y = element_text(size = 7), 
            axis.text.x = element_text(size = 12), 
            legend.title = element_text(size = 12), 
            legend.text = element_text(size = 12) 
        )
    return(p)
}

run_go <- function(genes_go, pvalue=0.05) {
    ego <- enrichGO(
        gene = genes_go,
        OrgDb = "org.Hs.eg.db",
        ont = "BP",
        pAdjustMethod = "BH",
        keyType = "SYMBOL",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        readable = TRUE
    )
    if (!is.null(ego) && "qvalue" %in% colnames(as.data.frame(ego))) {                                    
        ego_df <- as.data.frame(ego)                                    

        sorted_ego <- ego_df[order(ego_df$qvalue, decreasing = TRUE), ]   
        sorted_ego <- sorted_ego[sorted_ego$pval <= pvalue, ]
        sorted_ego <- sorted_ego[!is.na(sorted_ego$pval), ]  
        return(sorted_ego)                                                 
    } else {
        stop("Error: The 'ego' object is NULL or does not contain a 'qvalue' column.")      
    }
}


plot_go <- function(sorted_ego, pvalue=0.05) {
  plot_data <- sorted_ego %>%
    filter(!is.na(pvalue), pvalue <= pvalue) %>% 
    arrange(desc(Count)) %>%
    slice_head(n = 20) %>%
    mutate(
      GO_term         = fct_reorder(Description, Count),
      minus_log10_pval = -log10(pvalue)
    )
  
  p <- ggplot(plot_data, aes(x = minus_log10_pval, y = GO_term)) +
    geom_point(aes(size = Count), color = "steelblue") +
    scale_size_continuous(name = "Gene count", range = c(3, 10)) + 
    labs(
      title   = "GO Analysis",
      x       = expression(-log[10](p-value)),
      y       = NULL,
      caption = NULL
    ) +
    theme_minimal(base_size = 14) + 
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold"), 
      axis.text.y  = element_text(size = 16), 
      axis.text.x  = element_text(size = 16), 
      legend.title = element_text(size = 12), 
      legend.text  = element_text(size = 12) 
    )
  
  return(p)
}



load_msigdbr_data <- function() {

    filtered_genesets <- msigdbr(species = "Homo sapiens", category = "C2") %>%
        filter(grepl("^CP", gs_subcat)) %>%
        group_by(gs_name) %>%
        filter(n() >= 5 & n() <= 500) %>%
        ungroup()

    filtered_genesets$gene_symbol <- toupper(filtered_genesets$gene_symbol)
    m_list <- filtered_genesets %>% split(x = filtered_genesets$gene_symbol, f = filtered_genesets$gs_name)
    return(m_list)
}

get_lineage_plots <- function(OUTPUTS, condition_name, assocRes, genes, lineage_name, m_list, gsea_pvalue=0.05, go_pvalue=0.05) {
    OUTPUTS[["Lineage_Plots"]][[lineage_name]] <- list("GSEA" = list(), "GO" = list())
    OUTPUTS[["Lineage_Plots"]][[lineage_name]][["GSEA"]] <- list("Plots" = list(), "Data" = list())
    OUTPUTS[["Lineage_Plots"]][[lineage_name]][["GO"]] <- list("Plots" = list(), "Data" = list())

    message(paste0("Applying GSEA and GO analysis for lineage ", lineage_name))

    if (condition_name=="all") {
        wald_test_column <- paste0("waldStat")
    } else{
        wald_test_column <- paste0("waldStat_", lineage_name, "_condition", condition_name)
    }

    gsea_res <- run_gsea(assocRes, genes, wald_test_column, m_list)
    fgsea_res <- gsea_res$fgseaRes
    OUTPUTS[["Lineage_Plots"]][[lineage_name]][["GSEA"]][["Data"]] <- fgsea_res
    outfile_fgsea <- paste0("./out/all_samples_fgsea.xlsx")
    write.xlsx(fgsea_res, outfile_fgsea, rowNames = TRUE)
    gsea_plot <- plot_fgsea(fgsea_res, pvalue=gsea_pvalue)
    OUTPUTS[["Lineage_Plots"]][[lineage_name]][["GSEA"]][["Plots"]] <- gsea_plot
            
    sorted_ego <- run_go(genes)
    OUTPUTS[["Lineage_Plots"]][[lineage_name]][["GO"]][["Data"]] <- sorted_ego
    outfile_go <- paste0("./out/all_samples_GO.xlsx")
    write.xlsx(sorted_ego, outfile_go, rowNames = TRUE)
    go_plot <- plot_go(sorted_ego, pvalue=go_pvalue)
    OUTPUTS[["Lineage_Plots"]][[lineage_name]][["GO"]][["Plots"]] <- go_plot

    return(OUTPUTS)
}


get_cluster_plots <- function(OUTPUTS, assocRes, condition_name, lineage_name, cluster_id, m_list, lineage_clusters, gsea_pvalue, go_pvalue) {
    
    current_cluster_genes <- lineage_clusters[[cluster_id]]  
     

    if (length(current_cluster_genes) < 10) {
        message(paste0("Skipping condition=", condition_name, 
                                ", lineage=", lineage_name, ", 
                                cluster=", cluster_id, 
                                " since the number of genes within this cluter is < 10, 
                                which is too low for GSEA and GO analysis"))
        return(OUTPUTS)
    }

    OUTPUTS[[lineage_name]][[cluster_id]] <- list("GSEA" = list(), "GO" = list(), "Genes" = c())
    OUTPUTS[[lineage_name]][[cluster_id]][["Genes"]] <- current_cluster_genes  
    OUTPUTS[[lineage_name]][[cluster_id]][["GSEA"]] <- list("Plots" = list(), "Data" = list())
    OUTPUTS[[lineage_name]][[cluster_id]][["GO"]] <- list("Plots" = list(), "Data" = list())

    message(paste0("Applying GSEA and GO analysis for ", condition_name, " lineage ", lineage_name, " cluster ", cluster_id))
    if (condition_name=="all") {
        wald_test_column <- paste0("waldStat")
    } else{
        wald_test_column <- paste0("waldStat_", lineage_name, "_condition", condition_name)
    }
    
    gsea_res <- run_gsea(assocRes, lineage_clusters[[cluster_id]], wald_test_column, m_list)
    fgsea_res <- gsea_res$fgseaRes
    OUTPUTS[[lineage_name]][[cluster_id]][["GSEA"]][["Data"]] <- fgsea_res
    
    outfile_fgsea <- paste0("./out/", condition_name, "_", lineage_name, "_", cluster_id, "_fgsea.xlsx")
    write.xlsx(fgsea_res, outfile_fgsea, rowNames = TRUE)
    gsea_plot <- plot_fgsea(fgsea_res, pvalue = gsea_pvalue)
    OUTPUTS[[lineage_name]][[cluster_id]][["GSEA"]][["Plots"]] <- gsea_plot
            
    sorted_ego <- run_go(lineage_clusters[[cluster_id]])
    OUTPUTS[[lineage_name]][[cluster_id]][["GO"]][["Data"]] <- sorted_ego
    outfile_go <- paste0("./out/", condition_name, "_", lineage_name, "_", cluster_id, "_GO.xlsx")
    write.xlsx(sorted_ego, outfile_go, rowNames = TRUE)
    go_plot <- plot_go(sorted_ego, pvalue = go_pvalue)
    OUTPUTS[[lineage_name]][[cluster_id]][["GO"]][["Plots"]] <- go_plot
    return(OUTPUTS)
}



######################################################################################

lineage_gene_collection <- function(OUTPUTS, assocRes, condition_name, LINEAGE_INDICES) {
    # OUTPUTS: A list with lineage names as names
    # assocRes: Association test result matrix
    # condition_name: Can be 2_percent, 6_percent, Ctrl, or all
    # Returns: 
    #    lineage_genes: a list - contains lineage name as names and significant genes obtained
    #                            from association test results using corresponding p values columns
    #    genes: genes vector - contains combination of significant genes obtained 
    #                          from each lineage p-value column from association test
    #          

    message(paste0("Pulling genes associated with lineages for ", condition_name))
    genes_res <- get_lineage_genes(assocRes, condition_name, LINEAGE_INDICES)
    lineage_genes <- list()
    for (lineage_name in names(LINEAGE_INDICES)) {
        OUTPUTS[[lineage_name]] <- list()
        lineage_genes[[lineage_name]] <- genes_res[[lineage_name]]
    }

    genes <- Reduce(union, lineage_genes)
    return(list(
            "genes" = genes,
            "lineage_genes" = lineage_genes
        )
    )
}

generate_yhatsmooth <- function(sce, genes, condition_name) {
    # sce: fitgam output single cell experiment object
    # genes: genes vector - contains combination of significant genes obtained 
    #                       from each lineage p-value column from association test
    # condition_name: Can be 2_percent, 6_percent, Ctrl, or all
    # Returns: 
    #    heatSmooth_res: a list - contains lineage names and yhatSmooth_condition_obj as names, contains following
    #         heatSmooth_res['yhatSmooth_condition_obj']: scaled yhatSmooth data columns for specific condition
    #         heatSmooth_res['lineage1'], heatSmooth_res['lineage2'].....: heatmap plots for each lineage obtained from get_heatSmooth

    message(paste0("Running predictSmooth for ", condition_name))
    yhatSmooth <- predictSmooth(sce, gene = genes, nPoints = 50, tidy = FALSE)
    if (condition_name == "all") {
        keep_cols <- colnames(yhatSmooth)
    } else {
        keep_cols <- grepl(paste0("condition", condition_name), colnames(yhatSmooth))
    }
    
    yhatSmooth_condition_obj <- yhatSmooth[, keep_cols, drop = FALSE]
    #yhatSmoothSdataMatrix <- t(scale(t(yhatSmooth_condition_obj)))
    heatSmooth_res <- get_heatSmooth(yhatSmooth_condition_obj, LINEAGE_INDICES)
    heatSmooth_res[["yhatSmooth_dataMatrix"]] <- yhatSmooth_condition_obj
    return(heatSmooth_res)
}

generate_heatmap <- function(yhatSmooth_obj) {
    # yhatSmooth_obj: heatmap plots for each lineage obtained from get_heatSmooth
    # Returns: 
    #    heatmap_grid_plot: combined heatmap of all lineages

    heatmap_grid_plot <- arrangeGrob(yhatSmooth_obj[["lineage1"]][[4]], 
                                    yhatSmooth_obj[["lineage2"]][[4]], 
                                    yhatSmooth_obj[["lineage3"]][[4]], 
                                    yhatSmooth_obj[["lineage4"]][[4]], 
                                    yhatSmooth_obj[["lineage5"]][[4]], 
                                    ncol = 5)
    return(heatmap_grid_plot)
}


get_cluster_annotation_and_extract_genes <- function(heatSmooth_list,
                                                    LINEAGE_INDICES,
                                                    clustering_range
                                                    ) {
    # heatSmooth_list: output of generate_yhatsmooth which contains,
    #         heatSmooth_list['yhatSmooth_dataMatrix']: yhatSmooth data columns of the specific sample condition under test
    #         heatSmooth_list['lineage1'], heatSmooth_res['lineage2'].....: heatmap plots for each lineage obtained from get_heatSmooth
    # Returns: 
    #     optimal_k_list: optimal number of cluster determined by applying kmeans and calculating silhoutte score for cluster counts.
    #                     cluster count with highest silhoutte score is selected as optimal cluster count.
    #     cluster_genes: list of list with lineage name as names. Here, each sublist is a list with cluster id as name and genes as vector.
    #      
    cluster_genes <- list("lineage1" = list(), 
                          "lineage2" = list(),
                          "lineage3" = list(),
                          "lineage4" = list(),
                          "lineage5" = list())

    optimal_k_list <- list()
    lineage_processing_ret_list <- list()
    for (lineage_name in names(LINEAGE_INDICES)) {
        lineage_processing_ret_list[[lineage_name]] <- get_lineage_cluster(cluster_genes, 
                                                                            heatSmooth_list[[lineage_name]], 
                                                                            heatSmooth_list[["yhatSmooth_dataMatrix"]], 
                                                                            lineage_label=lineage_name, 
                                                                            LINEAGE_INDICES, 
                                                                            clustering_range)

        cluster_genes <- lineage_processing_ret_list[[lineage_name]][["cluster_genes"]]
        optimal_k_list[[lineage_name]] <- lineage_processing_ret_list[[lineage_name]][["optimal_k"]]
    }
    
    return(list("cluster_genes" = cluster_genes,
                "optimal_k_list" = optimal_k_list
        )
    )
}

get_cluster_annotation <- function(heatSmooth_list, optimal_k_list, LINEAGE_INDICES) {
    # heatSmooth_list: list of heatmaps (pheatmaps) with clusters for all lineages
    # optimal_k_list: optimal number of clusters determined by kmeans slihoutte score
    # LINEAGE_INDICES: list of lineage column indices in output obtained from predictSmooth function. Contains lineange names as names.
    # Return:`
    #    ann_combined: Gene clusters, with their corresponding yhatSmooth_rets[['yhatSmooth_dataMatrix']] values, for each lineage 

    ann_lineage_list <- list()
    for (lineage_name in names(LINEAGE_INDICES)) {
        cl_lineage_obj <- sort(cutree(heatSmooth_list[[lineage_name]]$tree_row, k = optimal_k_list[[lineage_name]]))
        ann_lineage_list[[lineage_name]] <- data.frame(cl_lineage_obj)
    }
    ann_combined <- do.call(cbind, ann_lineage_list)
    colnames(ann_combined) <- names(ann_lineage_list)
    return(ann_combined)
}

association_test_workflow <- function(sce, assocRes, condition_name, LINEAGE_INDICES, clustering_range, gsea_pvalue=0.05, go_pvalue=0.05) {
    # sce: fitgam single cell object
    # assocRes : association test results
    # condition_name: experiment condition under test - example: "2_percent", "6_percent", "Ctrl".
    #                   Use condition_name="all" when running it on fitgam results obtained from its global test (no conditional covariates used) 
    # LINEAGE_INDICES: Column indices of each lineage in predictSmooth results. 
    #                   Workflow uses nPoints=50 at the moment, which would lead to 50 columns for each indices
    #                   When suplying fitgam sce object obtained without using covariates, each lineage produces 50 columns, although
    #                   when covariates are used (condition_name="2_percent", "6_percent", or "Ctrl"), each lineage will have multiple covariates
    #                   Leading to 50 columns for each lineage-condition pairs.
    # clustering_range: Range of cluster counts that needs to be evluated to obtain optimal number of clusters. Values between range 2 to 10 works can be picked here.
    # gsea_pvalue: p value used to pick significantly enriched pathways in GSEA analysis 
    # go_pvalue: p value used to pick significantly enriched GO term in GO analysis
    # Returns: example:
    #       OUTPUTS: list(
    #                       "Heatmap"= ...,
    #                       "Lineage_Plots"=(
    #                                        "lineage1"=(
    #                                                    "GSEA"=...,
    #                                                    "GO"=...,
    #                                                    ),
    #                                        "lineage2":..., 
    #                                        "lineage3":...,  
    #                                        "lineage4":...,
    #                                        "lineage5":...,
    #                                        ),
    #                       "lineage1"=(
    #                                    "GSEA"=...,
    #                                    "GO"=...,
    #                                   ),
    #                        "lineage2":..., 
    #                        "lineage3":...,  
    #                        "lineage4":...,
    #                        "lineage5":...,
    #                        ),
    #                   )
    #                   
    
    OUTPUTS <- list()
    gene_results <- lineage_gene_collection(OUTPUTS, assocRes, condition_name, LINEAGE_INDICES)
    lineage_genes = gene_results$lineage_genes
    genes = gene_results$genes
    yhatSmooth_rets = generate_yhatsmooth(sce, genes, condition_name)
    heatmap_grid_plot <- generate_heatmap(yhatSmooth_rets)

    message(paste0("Clustering lineages for ", condition_name))
    workflow_ret <- get_cluster_annotation_and_extract_genes(yhatSmooth_rets, LINEAGE_INDICES, clustering_range)

    cluster_genes <- workflow_ret[["cluster_genes"]]
    optimal_k_list <- workflow_ret[["optimal_k_list"]]
    heatSmooth_res_with_clusters <- add_optimal_cluster_to_heatSmooth(yhatSmooth_rets[['yhatSmooth_dataMatrix']], 
                                                                    optimal_k_list, 
                                                                    LINEAGE_INDICES)

    heatmap_grid_plot_with_clusters <- generate_heatmap(heatSmooth_res_with_clusters)
    OUTPUTS[["Heatmap"]] <- heatmap_grid_plot_with_clusters

    ann_combined <- get_cluster_annotation(heatSmooth_res_with_clusters, optimal_k_list, LINEAGE_INDICES)
    outfile <- paste0("./out/", condition_name, "_association_test_cluster_annotations.xlsx")
    write.xlsx(ann_combined, outfile, rowNames = TRUE)

    message(paste0("Applying GSEA and GO analysis for ", condition_name))
    m_list <- load_msigdbr_data()
    OUTPUTS[["Lineage_Plots"]] <- list()
    for (lineage_name in names(cluster_genes)) {
        lineage_genes_vector <- lineage_genes[[lineage_name]]
        OUTPUTS <- get_lineage_plots(OUTPUTS, condition_name, assocRes, lineage_genes_vector, lineage_name, m_list, gsea_pvalue, go_pvalue)
        lineage_clusters <- cluster_genes[[lineage_name]]
        for (cluster_id in names(lineage_clusters)) {
            OUTPUTS <- get_cluster_plots(OUTPUTS, assocRes, condition_name, lineage_name, cluster_id,  m_list, lineage_clusters, gsea_pvalue, go_pvalue)
        }
    }
    return(OUTPUTS)
}


collect_fgsea_results <- function(OUTPUTS, condition_name) {
  results_list <- list()
  for (lineage_name in names(OUTPUTS)) {
    if (lineage_name == "Heatmap") next
    if (lineage_name == "Lineage_Plots") next
    clusters_obj <- OUTPUTS[[lineage_name]]
    
    for (cluster_id in names(clusters_obj)) {
        fgsea_data <- clusters_obj[[cluster_id]][["GSEA"]][["Data"]]
        if (is.null(fgsea_data) || !is.data.frame(fgsea_data)) next

        if (dim(fgsea_data)[1]==0) next

        df <- as.data.frame(fgsea_data)
        df$condition <- condition_name
        df$lineage   <- lineage_name
        df$cluster   <- cluster_id
        results_list[[length(results_list) + 1]] <- df
        }
    }
  
    all_results <- dplyr::bind_rows(results_list)
    return(all_results)
}
