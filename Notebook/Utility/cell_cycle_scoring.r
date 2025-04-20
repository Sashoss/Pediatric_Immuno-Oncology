suppressPackageStartupMessages({
    library(Seurat)
    library(tools) 
    library(AnnotationDbi)
    library(org.Hs.eg.db)  
    library(dplyr)
    library(AnnotationHub)
})


cell_cycle_scoring <- function(seurat_obj, cc_file="./Data/CC_Homo_sapiens.csv") {
    cell_cycle_genes <- read.csv(file = cc_file)
    ah <- AnnotationHub()
    ahDb <- query(ah,  pattern = c("Homo sapiens", "EnsDb"), ignore.case = TRUE)
    id <- ahDb %>% mcols() %>% rownames() %>% tail(n = 1)
    edb <- ah[[id]]
    annotations <- genes(edb, return.type = "data.frame")
    annotations <- annotations %>%
        dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

    cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))
    s_genes <- cell_cycle_markers %>%
        dplyr::filter(phase == "S") %>%
        pull("gene_name")
        
    g2m_genes <- cell_cycle_markers %>%
        dplyr::filter(phase == "G2/M") %>%
        pull("gene_name")

    seurat_obj <- CellCycleScoring(seurat_obj, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)
    seurat_obj$CC.Difference <- seurat_obj$S.Score - seurat_obj$G2M.Score

    return(seurat_obj)
}

