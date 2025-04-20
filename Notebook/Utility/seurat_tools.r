
clean_seurat <- function(seuratObj) {
    seuratObj <- DietSeurat(seuratObj, counts = TRUE, data = TRUE, scale.data = FALSE, 
                            dimreducs = NULL, graphs = FALSE)
    seuratObj@assays[["RNA"]] <- seuratObj@assays[["originalexp"]]
    seuratObj@assays[["originalexp"]] <- NULL
    DefaultAssay(seuratObj) <- "RNA"
    seuratObj[["nFeature_RNA"]] <- apply(seuratObj@assays$RNA@counts, 2, function(x) sum(x > 0))
    return(seuratObj)
}


get_module_score <- function(seurat_obj, marker_list) {
    score_names <- c()
    for (cell_type in names(marker_list)) {
        genes <- marker_list[[cell_type]]
        if (length(genes) <= 3) {
            next
        }
        valid_genes <- genes[genes %in% rownames(seurat_obj)]

        if (length(valid_genes) >= 2) {
            score_prefix <- cell_type
            
            seurat_obj <- AddModuleScore(
                                    seurat_obj,
                                    features = list(valid_genes),
                                    name = score_prefix
                        )

            added_col <- grep(paste0("^", score_prefix, "[0-9]+$"), colnames(seurat_obj[[]]), value = TRUE)
            score_names <- c(score_names, added_col)
        } else {
            warning(paste("Skipping", cell_type, "- insufficient valid genes found in Seurat object"))
        }
    }

    return(list(
            "seurat_obj"=seurat_obj,
            "scores_names"=score_names)
            )
}


