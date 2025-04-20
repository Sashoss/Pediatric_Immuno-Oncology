


suppressPackageStartupMessages({
    library(Seurat)
    library(DoubletFinder)
    library(dplyr)
    library(optparse)
})



callDoubletFinder <- function(seuratObj, min_pc, doublet_rate=0.075, pN=0.25, res=0.05, n_cores=1) {

    seuratObj <- NormalizeData(seuratObj)
    seuratObj <- FindVariableFeatures(seuratObj)
    seuratObj <- ScaleData(seuratObj)
    seuratObj <- RunPCA(seuratObj, npcs = min_pc, features = VariableFeatures(seuratObj))

    seuratObj <- RunUMAP(seuratObj, dims = 1:min_pc)
    seuratObj <- FindNeighbors(seuratObj, dims = 1:min_pc)              
    seuratObj <- FindClusters(seuratObj, resolution = res)

    sweep.list <- paramSweep(seuratObj, PCs = 1:min_pc, num.cores = n_cores) 
    sweep.stats <- summarizeSweep(sweep.list)
    bcmvn <- find.pK(sweep.stats)

    bcmvn_max <- bcmvn[which.max(bcmvn$BCmetric), ]
    optimal_pk <- as.numeric(as.character(bcmvn_max$pK))

    annotations <- seuratObj@meta.data$seurat_clusters
    homotypic_prop <- modelHomotypic(annotations) 
    nExp_poi <- round(doublet_rate * nrow(seuratObj@meta.data)) 
    nExp_poi_adj <- max(0, round(nExp_poi * (1 - homotypic_prop)))

    seuratObj <- doubletFinder(seu = seuratObj, 
                            PCs = 1:min_pc, 
                            pN = pN, 
                            pK = optimal_pk,
                            nExp = nExp_poi_adj)

    seuratObj@meta.data <- seuratObj@meta.data %>%
                            rename_with(~ paste0("DFC_", seq_along(.)), matches("^DF.classifications")) %>%
                            rename_with(~ paste0("pANN_", seq_along(.)), matches("^pANN_"))

    return(seuratObj)
}
