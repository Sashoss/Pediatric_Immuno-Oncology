suppressPackageStartupMessages({
    library(Seurat)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(Matrix)
})


convert_ensembl_to_symbol <- function(ensembl_ids) {
    ensembl_ids_clean <- sub("\\..*", "", ensembl_ids)
    mapIds(org.Hs.eg.db,
        keys = ensembl_ids_clean,
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first")
}


update_gene_names <- function(seurat_obj) {
    genes <- rownames(seurat_obj)
    if (!any(grepl("^ENSG", genes))) {
        message("Gene names already appear to be symbols. No conversion needed.")
        return(seurat_obj)
    }

    message("Converting Ensembl IDs to gene symbols...")
    gene_symbols <- convert_ensembl_to_symbol(genes)

    valid <- !is.na(gene_symbols) & gene_symbols != ""
    if (sum(valid) == 0) {
        warning("No valid gene symbols found. Returning original object.")
        return(seurat_obj)
    }

    seurat_obj <- seurat_obj[valid, ]
    new_names <- gene_symbols[valid]
    if (any(duplicated(new_names))) {
        counts <- GetAssayData(seurat_obj, slot = "counts")
        counts_agg <- aggregate.Matrix(counts,
                                        groupings = new_names,
                                        fun = "sum")
        seurat_obj <- CreateSeuratObject(counts = counts_agg,
                                        meta.data = seurat_obj@meta.data,
                                        min.cells = 3,
                                        min.features = 200)
    } else {
        rownames(seurat_obj) <- new_names
    }

    seurat_obj[["nFeature_RNA"]] <- apply(GetAssayData(seurat_obj, slot = "counts"), 2, function(x) sum(x > 0))
    return(seurat_obj)
}


