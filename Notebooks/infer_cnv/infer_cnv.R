# set up

library(Seurat)
library(tidyverse)
library(infercnv)


# Merged data
seurat_rds_path <- "./out/integrated_seurat_clusters.rds"
seurat_obj <- readRDS(seurat_rds_path)
Idents(seurat_obj) <- seurat_obj$RNA_snn_res.0.05

# make count matrix from seurat object
matrix <- GetAssayData(seurat_obj,
                       layer = "counts",
                       assay = "RNA") %>%
            as.data.frame()


gene_order <- read.table("./in/inferCNV_inputs/gencode_v19_gene_pos.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
genes_in_counts <- rownames(matrix)
genes_in_order <- gene_order$V1

common_genes <- intersect(genes_in_counts, genes_in_order)
matrix <- matrix[common_genes, ]
write.table(matrix, 
            file = "./in/inferCNV_inputs/merged_gene_matrix.txt", 
            append = FALSE, 
            sep = "\t", 
            row.names = TRUE, 
            col.names = TRUE)

gene_order <- gene_order[gene_order$V1 %in% common_genes, ]
write.table(
    gene_order,
    file = "./in/inferCNV_inputs/inferCNV_gene_order.txt",
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
)



# make annotation with Seurat clusters

annotation <- Idents(merged) %>% as.matrix()

write.table(annotation, 

            file = "./in/inferCNV_inputs/merged_annotation.txt", 

            append = FALSE, 

            sep = "\t", 

            row.names = TRUE, 

            col.names = FALSE)



gc()


-----------------------------------------------------------------------------------------------------------------------------

run_infercnv.r

library(Seurat)
library(infercnv)
library(ggplot2)

options(scipen = 100) #to prevent scientific notation issues that can interfere with the hierarchical clustering process.


setwd("/home/gdwanglab/axk201/drissilab/Single_cell_Oct2024/Analysis/Dec_Workflow/cluster")

# create InferCNV object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="./in/inferCNV_inputs/merged_gene_matrix.txt",
                                    annotations_file="./in/inferCNV_inputs/merged_annotation.txt",
                                    delim="\t",
                                    gene_order_file="./in/inferCNV_inputs/inferCNV_gene_order.txt",
                                    ref_group_names=c("2"))

# run InferCNV
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="./out/inferCNV_output",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T,
                             num_threads = 24)

# Plot CNV Heatmap
p <- infercnv::plot_cnv(
    infercnv_obj,
    output_format = "png",
    title = "inferCNV Heatmap",
    cluster_by_groups = TRUE,
    plot_num_markers = 100,
    legend_title = "CNV Status",
    draw_lines = TRUE
)


ggsave(p, "./out/infercnv_heatmap.png")


