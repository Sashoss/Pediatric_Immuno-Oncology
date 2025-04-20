suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(infercnv)
    library(dplyr)
    library(ggplot2)
})



create_infercnv_input <- function(seurat_obj, 
                                gencode_path="./in/inferCNV_inputs/gencode_v19_gene_pos.txt", 
                                output_dir="./out/infercnv_input"
                                ) {

    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    matrix <- GetAssayData(seurat_obj,
                        layer = "counts",
                        assay = "RNA") %>%
                                    as.data.frame()

    gene_order <- read.table(file.path(gencode_path), header = FALSE, sep = "\t", stringsAsFactors = FALSE)

    genes_in_counts <- rownames(matrix)
    genes_in_order <- gene_order$V1
    common_genes <- intersect(genes_in_counts, genes_in_order)

    matrix <- matrix[common_genes, ]
    write.table(matrix, 
                file = file.path(output_dir, "merged_gene_matrix.txt"), 
                append = FALSE, 
                sep = "\t", 
                row.names = TRUE, 
                col.names = TRUE
    )

    gene_order <- gene_order[gene_order$V1 %in% common_genes, ]
    write.table(
                gene_order,
                file = file.path(output_dir, "inferCNV_gene_order.txt"),
                sep = "\t",
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE
    )

    annotation <- Idents(seurat_obj) %>% as.matrix()
    write.table(annotation, 
            file = file.path(output_dir, "merged_annotation.txt"), 
            append = FALSE, 
            sep = "\t", 
            row.names = TRUE, 
            col.names = FALSE)

}



run_infercnv <- function(
                    raw_counts_matrix="./in/inferCNV_inputs/merged_gene_matrix.txt",
                    annotations_file="./in/inferCNV_inputs/merged_annotation.txt",
                    gene_order_file="./in/inferCNV_inputs/inferCNV_gene_order.txt",
                    ref_group_names==c("NK cells"),
                    cutoff=0.1,
                    num_threads=24,
                    output_dir="./out/inferCNV_output"
                ) {

    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=raw_counts_matrix,
                                    annotations_file=annotations_file,
                                    delim="\t",
                                    gene_order_file=gene_order_file,
                                    ref_group_names=ref_group_names
                    )

    infercnv_obj = infercnv::run(infercnv_obj,
                            cutoff=cutoff,  
                            out_dir=output_dir,  
                            cluster_by_groups=T,   
                            denoise=T,
                            HMM=T,
                            num_threads=num_threads)

    return(infercnv_obj)
}